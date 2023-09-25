# -*- coding: utf-8 -*-
"""
Created August 2021

@author: Poornima Sivanandam

Script to process DJI Zenmuse P1 (gimbal 1) and MicaSense RedEdge-MX/Dual (gimbal 2) images captured simultaneously
using the Matrice 300 RTK drone system.

Assumption that folder structure is as per the TERN protocols:
Data |	Path | Example
Raw data |	<plot>/YYYYMMDD/imagery/<sensor>/level0_raw/ |	SASMDD0001/20220519/imagery/rgb/level0_raw
Data products |	<plot>/YYYYMMDD/imagery/<sensor>/level1_proc/	| SASMDD0001/20220519/imagery/multispec/level1_proc
Metashape project |	plot/YYYYMMDD/imagery/metashape| SASRIV0001/20220516/imagery/metashape/
DRTK logs | plot/YYYYMMDD/drtk/

Raw data paths can be overriden using 'Optional Inputs'.

Required Input:
    -crs "<EPSG code for target projected coordinate reference system. Also used in MicaSense position interpolation>"
    Example: -crs "7855"
    See https://epsg.org/home.html

Optional Inputs:
    1. -multispec "path to multispectral level0_raw folder containing raw data"
        Default is relative to project location: ../multispec/level0_raw/
    2. -rgb "path to RGB level0_raw folder which also has the MRK file(s)"
        Default is relative to project location: ../rgb/level0_raw/
    3. -smooth "<low/medium/high>"
        Strength value to smooth RGB model. Default is low.
        Low: for low-lying vegetation (grasslands, shrublands), Medium and high: as appropriate for forested sites.
    4. When P1 (RGB camera) coordinates have to be blockshifted:
        - Path to file containing DRTK init and AUSPOS cartesian coords passed using "-drtk <path to file>".

Summary:
    * Add RGB and multispectral images.
    * Stop script for user input on calibration images.
    * When 'Resume Processing' is clicked complete the processing workflow.

"""

import argparse
import math
import collections
import numpy as np
import Metashape
import os
import sys
import exifread
from collections import defaultdict
from upd_micasense_pos import ret_micasense_pos
import importlib
import upd_micasense_pos

importlib.reload(upd_micasense_pos)
from pathlib import Path

# Note: External modules imported were installed through:
# "C:\Program Files\Agisoft\Metashape Pro\python\python.exe" -m pip install <modulename>
# See M300 data processing protocol for more information.

# Metashape Python API updates in v2.0
METASHAPE_V2_PLUS = False
found_version = Metashape.app.version.split('.')  # e.g. 2.0.1
if int(found_version[0]) >= 2:
    METASHAPE_V2_PLUS = True

###############################################################################
# Constants
###############################################################################
GEOG_COORD = collections.namedtuple('Geog_CS', ['lat_decdeg', 'lon_decdeg', 'elliph'])

SOURCE_CRS = Metashape.CoordinateSystem("EPSG::4326")  # WGS84

CONST_a = 6378137  # Semi major axis
CONST_inv_f = 298.257223563  # Inverse flattening 1/f WGS84 ellipsoid
# Chunks in Metashape
CHUNK_RGB = "rgb"
CHUNK_MULTISPEC = "multispec"

IMG_QUAL_THRESHOLD = 0.7

DICT_SMOOTH_STRENGTH = {'low': 50, 'medium': 100, 'high': 200}

# Lever-arm offsets for different sensors on *Matrice 300*
# TODO: update this for other sensors and drone platforms
P1_GIMBAL1_OFFSET = (0.087, 0.0, 0.0)

# Measure lever-arm offsets (X, Y, Z) from the single gimbal position to the ‘master’ camera (by default, the lowest wavelength)
# In Metashape, the offsets are positive with respect to the actual camera positions. 
# See Metashape manual or TERN RGB Multispectral processing protocol for details.
offset_dict = defaultdict(dict)
offset_dict['RedEdge-M']['Red'] = (-0.097, -0.03, -0.06)
offset_dict['RedEdge-M']['Dual'] = (-0.097, 0.02, -0.08)
offset_dict['RedEdge-P']['Red'] = (0,0,0)
offset_dict['RedEdge-P']['Dual'] = (0,0,0)

###############################################################################
# Function definitions
###############################################################################
def cartesian_to_geog(X, Y, Z):
    """
    Author: Poornima Sivanandam
    Convert Cartesian coordinates to geographic coordinates using WGS84 ellipsoid.
    Return Lat, Lon, ellipsoidal height as a named tuple.
    Calculations from Transformation_Conversion.xlsx at https://github.com/icsm-au/DatumSpreadsheets
    """
    f = 1 / CONST_inv_f
    e_sq = 2 * f - f ** 2
    p = math.sqrt(X ** 2 + Y ** 2)
    r = math.sqrt(p ** 2 + Z ** 2)
    mu = math.atan((Z / p) * (1 - f) + (e_sq * CONST_a) / r)

    lat_top_line = Z * (1 - f) + e_sq * CONST_a * math.sin(mu) ** 3
    lat_bottom_line = (1 - f) * (p - e_sq * CONST_a * math.cos(mu) ** 3)

    lon = math.atan(Y / X)
    lat = math.atan(lat_top_line / lat_bottom_line)

    if (lon < 0):
        tmp_lon = lon + math.pi
    else:
        tmp_lon = lon

    lon_dec_deg = (tmp_lon / math.pi) * 180
    lat_dec_deg = (lat / math.pi) * 180

    ellip_h = p * math.cos(lat) + Z * math.sin(lat) - CONST_a * math.sqrt(1 - e_sq * math.sin(lat) ** 2)

    conv_coord = GEOG_COORD(lat_dec_deg, lon_dec_deg, ellip_h)

    return conv_coord


def find_files(folder, types):
    photo_list = list()
    for dir, subdir, file in os.walk(folder):
        for filename in file:
            if (filename.lower().endswith(types)):
                photo_list.append(os.path.join(dir, filename))
    return (photo_list)


def proc_rgb():
    """
    Author: Poornima Sivanandam
    Arguments: None
    Return: None
    Create: RGB orthomosaic in rgb/level1_proc or in Metashape project folder
        smoothed 3D model file in Metashape project folder
    Summary:
        * blockshift (optional through args)
        * convert to target CRS
        * Image Quality check
        * Apply GPS/INS offset for gimbal 1
        * Update Camera Accuracy settings for M300 RTK GNSS accuracy
        * Align images
        * Build dense cloud
        * Build model, decimate and smooth (use args)
        * Export model (for multispec chunk)
        * Build and export orthomosaic
    """
    # If P1 positions are to be blockshifted, do the following:
    # - Read the .txt file and convert Cartesian coordinates to WGS84 Lat/Lon
    # - Calculate the difference and apply the shift directly to the cameras (Lon/Lat/Ellipsoidal height) in 'rgb' chunk
    # Convert coordinate system for Lat/Lon to target projected coordinate system

    chunk = doc.findChunk(dict_chunks[CHUNK_RGB])
    proj_file = doc.path
    blockshift_p1 = False

    if args.drtk is not None:
        blockshift_p1 = True
        DRTK_TXT_FILE = args.drtk
        print("P1 blockshift set")

        # read from txt/csv cartesian for RTK initial (line 1) and AUSPOS coords (line 2)
        with open(DRTK_TXT_FILE, 'r') as file:
            line = file.readline()
            split_line = line.split(',')
            drtk_field = cartesian_to_geog(float(split_line[0]), float(split_line[1]), float(split_line[2]))
            line = file.readline()
            split_line = line.split(',')
            drtk_auspos = cartesian_to_geog(float(split_line[0]), float(split_line[1]), float(split_line[2]))

        # calc difference
        diff_lat = round((drtk_auspos.lat_decdeg - drtk_field.lat_decdeg), 6)
        diff_lon = round((drtk_auspos.lon_decdeg - drtk_field.lon_decdeg), 6)
        diff_elliph = round((drtk_auspos.elliph - drtk_field.elliph), 6)
        P1_shift = Metashape.Vector((diff_lon, diff_lat, diff_elliph))

        print("Shifting P1 cameras by: " + str(P1_shift))

        # shift coordinates in the chunk
        for camera in chunk.cameras:
            if not camera.label == camera.master.label:
                continue
            if not camera.reference.location:
                continue
            else:
                camera.reference.location = camera.reference.location + P1_shift

    # Convert to projected coordinate system
    target_crs = Metashape.CoordinateSystem("EPSG::" + args.crs)
    for camera in chunk.cameras:
        if not camera.reference.location:
            continue
        camera.reference.location = Metashape.CoordinateSystem.transform(camera.reference.location, SOURCE_CRS,
                                                                         target_crs)

    chunk.crs = target_crs

    global P1_shift_vec
    if blockshift_p1:
        # Export updated positions as csv for debug purposes. Not used in script.
        chunk.exportReference(path=str(P1_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                              delimiter=",", items=Metashape.ReferenceItemsCameras)

        # If P1  blockshifted, pass vector for x, y, z shift of micasense image position
        P1_shift_vec = np.array([diff_lat, diff_lon, diff_elliph])
    else:
        P1_shift_vec = np.array([0.0, 0.0, 0.0])

    doc.save()

    #
    # Estimate image quality and remove cameras with quality < threshold
    #
    if METASHAPE_V2_PLUS:
        chunk.analyzeImages()
    else:
        chunk.analyzePhotos()
    low_img_qual = []
    low_img_qual = [camera for camera in chunk.cameras if (float(camera.meta["Image/Quality"]) < IMG_QUAL_THRESHOLD)]
    if low_img_qual:
        print("Removing cameras with Image Quality < %.1f" % IMG_QUAL_THRESHOLD)
        chunk.remove(low_img_qual)
    doc.save()

    #
    # GPS/INS offset
    #
    print(chunk.sensors[0].antenna.location_ref)
    print("Update GPS/INS offset for P1")
    chunk.sensors[0].antenna.location_ref = Metashape.Vector(P1_GIMBAL1_OFFSET)
    print(chunk.sensors[0].antenna.location_ref)

    #
    # Align Photos
    #
    print("Aligning Cameras")
    # change camera position accuracy to 0.1 m
    chunk.camera_location_accuracy = Metashape.Vector((0.10, 0.10, 0.10))

    # Downscale values per https://www.agisoft.com/forum/index.php?topic=11697.0
    # Downscale: highest, high, medium, low, lowest: 0, 1, 2, 4, 8
    # Quality:  High, Reference Preselection: Source
    chunk.matchPhotos(downscale=1, generic_preselection=False, reference_preselection=True,
                      reference_preselection_mode=Metashape.ReferencePreselectionSource)
    chunk.alignCameras()
    doc.save()

    #
    # Optimise Cameras
    #
    print("Optimise alignment")
    chunk.optimizeCameras()
    doc.save()

    #
    # Build Dense Cloud
    #
    # check if exists and reuse depthmap? # reuse_depth=True below
    # downscale: ultra, high, medium, low, lowest: 1, 2, 4, 8, 16
    print("Build dense cloud")
    # Medium quality. And default: mild filtering.
    chunk.buildDepthMaps(downscale=4)
    doc.save()

    if METASHAPE_V2_PLUS:
        chunk.buildPointCloud()
    else:
        chunk.buildDenseCloud()
    doc.save()

    #
    # Build Mesh
    #
    print("Build mesh")
    if METASHAPE_V2_PLUS:
        chunk.buildModel(surface_type=Metashape.HeightField, source_data=Metashape.PointCloudData,
                         face_count=Metashape.MediumFaceCount)
    else:
        chunk.buildModel(surface_type=Metashape.HeightField, source_data=Metashape.DenseCloudData,
                         face_count=Metashape.MediumFaceCount)
    doc.save()

    # Decimate and smooth mesh to use as orthorectification surface
    # Halve face count?
    chunk.decimateModel(face_count=len(chunk.model.faces) / 2)
    # Smooth model
    smooth_val = DICT_SMOOTH_STRENGTH[args.smooth]
    chunk.smoothModel(smooth_val)
    # Export model for use in micasense chunk
    model_file = Path(proj_file).parent / (Path(proj_file).stem + "_rgb_smooth_" + str(smooth_val) + ".obj")
    chunk.exportModel(path=str(model_file), crs=target_crs, format=Metashape.ModelFormatOBJ)

    #
    # Build and export orthomosaic
    #
    print("Build orthomosaic")
    chunk.buildOrthomosaic(surface_data=Metashape.DataSource.ModelData, refine_seamlines=True)
    doc.save()

    if chunk.orthomosaic:
        # Round resolution to 2 decimal places
        res_xy = round(chunk.orthomosaic.resolution, 2)

        # if rgb/ folder does not exist in MRK_PATH save orthomosaic in the project directory
        # else save ortho in rgb/level1_proc/
        p1_idx = MRK_PATH.find("rgb")
        if p1_idx == -1:
            dir_path = Path(proj_file).parent
            print("Cannot find rgb/ folder. Saving ortho in " + str(dir_path))
        else:
            # create p1/level1_proc folder if it does not exist
            dir_path = Path(MRK_PATH[:p1_idx + len("rgb")]) / "level1_proc"
            dir_path.mkdir(parents=True, exist_ok=True)

        # file naming format: <projname>_rgb_ortho_<res_in_m>.tif
        ortho_file = dir_path / (
                Path(proj_file).stem + "_rgb_ortho_" + str(res_xy).split('.')[1] + ".tif")

        compression = Metashape.ImageCompression()
        compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW  # default on Metashape
        compression.tiff_big = True
        compression.tiff_tiled = True
        compression.tiff_overviews = True

        chunk.exportRaster(path=str(ortho_file), resolution_x=res_xy, resolution_y=res_xy,
                           image_format=Metashape.ImageFormatTIFF,
                           save_alpha=False, source_data=Metashape.OrthomosaicData, image_compression=compression)
        print("Exported orthomosaic " + str(ortho_file))

    print("RGB chunk processing complete!")


def proc_multispec():
    """
    Author: Poornima Sivanandam
    Arguments: None
    Return: None
    Create: Multispec orthomosaic in multispec/level1_proc or in Metashape project folder
    Summary:
        * Interpolate micasense image position using p1 pos and timestamp.
        * Remove images that triggered outside p1 capture times
        * Image Quality check
        * Apply GPS/INS offset for gimbal 2
        * Set primary channel to NIR
        * Update Camera Accuracy settings for M300 RTK GNSS accuracy
        * Set raster transform to export relative reflectance in orthomosaic
        * Calibrate reflectance using both sun senors and panels
        * Align images
        * Build dense cloud
        * Import RGB smoothed model (see proc_rgb)
        * Build and export orthomosaic with raster transformed values (relative reflectance)
    """

    chunk = doc.findChunk(dict_chunks[CHUNK_MULTISPEC])

    target_crs = Metashape.CoordinateSystem("EPSG::" + args.crs)

    # Get image suffix of master camera
    camera = chunk.cameras[0]
    cam_master = camera.master.label.split('_')

    # file naming assumption: IMG_xxxx_suffixNum
    img_suffix_master = cam_master[2]

    print("Interpolate Micasense position based on P1 with blockshift" + str(P1_shift_vec))

    # inputs: paths to MRK file for P1 position, Micasense image path, image suffix for master band images, target CRS
    # returns output csv file with interpolated micasense positions
    ret_micasense_pos(MRK_PATH, MICASENSE_PATH, img_suffix_master, args.crs,
                      str(MICASENSE_CAM_CSV), P1_shift_vec)

    # Load updated positions in the chunk
    chunk.importReference(str(MICASENSE_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                          delimiter=",", crs=target_crs, skip_rows=1,
                          items=Metashape.ReferenceItemsCameras)
    doc.save()

    # ret_micasense_pos wrote Altitude = 0 (last column) for MicaSense images that triggered when P1 did not.
    # Create a list of cameras with Altitude = 0
    del_camera_names = list()

    # Only look at altitude of master band images
    for camera in chunk.cameras:
        if not camera.label == camera.master.label:
            continue
        if not camera.reference.location:
            continue
        if camera.reference.location.z == 0:
            del_camera_names.append(camera.label)

    # Delete images outside of P1 capture times
    print("Deleting MicaSense images that triggered outside P1 capture times")
    for camera in chunk.cameras:
        # Only calibration images are in a group. The following line is necessary to avoid NoneType error on other images
        if camera.group is not None:
            if camera.group.label == 'Calibration images':
                continue
        if camera.label in del_camera_names:
            chunk.remove(camera)

    # save project
    doc.save()


    # Set primary channel
    #
    # Get index of NIR band. Micasense Dual: NIR is sensors[9], and in RedEdge-M sensors[4]
    if cam_model == 'RedEdge-M':
        set_primary = "NIR"
    elif cam_model == 'RedEdge-P':
        set_primary = 'Panchro'
    for s in chunk.sensors:
        if s.label.find(set_primary) != -1:
            print("Setting primary channel to " + s.label)
            chunk.primary_channel = s.layer_index
            break


    # GPS/INS offset for master sensor
    #
    print("Updating Micasense GPS offset")
    chunk.sensors[0].antenna.location_ref = Metashape.Vector(MS_GIMBAL2_OFFSET)

    #
    # Set Raster Transform to calculate reflectance
    #
    print("Updating Raster Transform for relative reflectance")
    raster_transform_formula = []
    num_bands = len(chunk.sensors)
    if cam_model == 'RedEdge-M':
        for band in range(1, num_bands + 1):
            raster_transform_formula.append("B" + str(band) + "/32768")
    elif cam_model == 'RedEdge-P':
        # Skip Panchromatic band in multispec ortho.
        # Panchro band: wavelength: 634.5 nm, Band 5 in RedEdge-P Dual and Band 3 in RedEdge-P.
        if num_bands >= 10:
            PANCHRO_BAND = 5
        else:
            PANCHRO_BAND = 3
        for band in range(1, num_bands+1):
            if band != PANCHRO_BAND:
                raster_transform_formula.append("B" + str(band) + "/32768")

    chunk.raster_transform.formula = raster_transform_formula
    chunk.raster_transform.calibrateRange()
    chunk.raster_transform.enabled = True
    doc.save()

    #
    # Estimate image quality and remove cameras with quality < threshold
    #
    if METASHAPE_V2_PLUS:
        chunk.analyzeImages()
    else:
        chunk.analyzePhotos()
    low_img_qual = []
    low_img_qual = [camera.master for camera in chunk.cameras if (float(camera.meta["Image/Quality"]) < 0.5)]
    if low_img_qual:
        print("Removing cameras with Image Quality < %.1f" % 0.5)
        chunk.remove(list(set(low_img_qual)))
    doc.save()
    #
    #
    # Calibrate Reflectance
    #
    chunk.calibrateReflectance(use_reflectance_panels=True, use_sun_sensor=True)

    #
    # Align Photos
    #
    # change camera position accuracy to 0.1 m
    chunk.camera_location_accuracy = Metashape.Vector((0.10, 0.10, 0.10))

    # Downscale values per https://www.agisoft.com/forum/index.php?topic=11697.0
    # Downscale: highest, high, medium, low, lowest: 0, 1, 2, 4, 8
    # Quality:  High, Reference Preselection: Source
    chunk.matchPhotos(downscale=1, generic_preselection=False, reference_preselection=True,
                      reference_preselection_mode=Metashape.ReferencePreselectionSource)
    doc.save()
    print("Aligning cameras")
    chunk.alignCameras()
    doc.save()

    #
    # Optimise Cameras
    #
    print("Optimise alignment")
    chunk.optimizeCameras()
    doc.save()

    #
    # Build and export orthomosaic
    #
    # Import P1 model for use in orthorectification
    smooth_val = DICT_SMOOTH_STRENGTH[args.smooth]
    model_file = Path(proj_file).parent / (Path(proj_file).stem + "_rgb_smooth_" + str(smooth_val) + ".obj")
    chunk.importModel(path=str(model_file), crs=target_crs, format=Metashape.ModelFormatOBJ)

    print("Build orthomosaic")
    chunk.buildOrthomosaic(surface_data=Metashape.DataSource.ModelData, refine_seamlines=True)
    doc.save()

    if chunk.orthomosaic:
        # Round resolution to 2 decimal places
        res_xy = round(chunk.orthomosaic.resolution, 2)

        # if multispec/ folder does not exist in MICASENSE_PATH save in project directory
        # else save ortho in multispec/level1_proc/
        micasense_idx = MICASENSE_PATH.find("multispec")
        if micasense_idx == -1:
            dir_path = Path(proj_file).parent
            print("Cannot find " + "multispec/ folder. Saving ortho in " + str(dir_path))
        else:
            # create multispec/level1_proc/ folder if it does not exist
            dir_path = Path(MICASENSE_PATH[:micasense_idx + len("multispec")]) / "level1_proc"
            dir_path.mkdir(parents=True, exist_ok=True)

        # file naming format: <projname>_multispec_ortho_<res_in_m>.tif
        ortho_file = dir_path / (
                Path(proj_file).stem + "_" + "multispec_ortho_" + str(res_xy).split('.')[1] + ".tif")

        compression = Metashape.ImageCompression()
        compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW  # default on Metashape
        compression.tiff_big = True
        compression.tiff_tiled = True
        compression.tiff_overviews = True

        chunk.exportRaster(path=str(ortho_file), resolution_x=res_xy, resolution_y=res_xy,
                           image_format=Metashape.ImageFormatTIFF,
                           raster_transform=Metashape.RasterTransformValue,
                           save_alpha=False, source_data=Metashape.OrthomosaicData, image_compression=compression)
        print("Exported orthomosaic: " + str(ortho_file))

    print("Multispec chunk processing complete!")


def resume_proc():
    # Process RGB chunk
    proc_rgb()
    # Process multispec chunk
    proc_multispec()
    print("End of script")


############################################
##  Main code
############################################
print("Script start")

# Parse arguments and initialise variables
parser = argparse.ArgumentParser(
    description='Update camera positions in P1 and/or MicaSense chunks in Metashape project')
parser.add_argument('-crs',
                    help='EPSG code for target projected CRS for micasense cameras. E.g: 7855 for GDA2020/MGA zone 55',
                    required=True)
parser.add_argument('-multispec', help='path to multispectral level0_raw folder with raw images')
parser.add_argument('-rgb', help='path to RGB level0_raw folder that also has the MRK files')
parser.add_argument('-smooth', help='Smoothing strength used to smooth RGB mesh low/med/high', default="low")
parser.add_argument('-drtk', help='If RGB coordinates to be blockshifted, file containing \
                                                  DRTK base station coordinates from field and AUSPOS')

global args
args = parser.parse_args()
global MRK_PATH, MICASENSE_PATH

# Metashape project
global doc
doc = Metashape.app.document
proj_file = doc.path

# if Metashape project has not been saved
if proj_file == '':
    if args.rgb:
        proj_file = str(Path(args.rgb).parents[0] / "metashape_project.psx")
        print("Metashape project saved as %s" % proj_file)
        doc.save(proj_file)

if args.rgb:
    MRK_PATH = args.rgb
else:
    # Default is relative to project location: ../rgb/level0_raw/
    MRK_PATH = Path(proj_file).parents[1] / "rgb/level0_raw"
    if not MRK_PATH.is_dir():
        sys.exit("%s directory does not exist. Check and input paths using -rgb " % str(MRK_PATH))
    else:
        MRK_PATH = str(MRK_PATH)

# TODO update when other sensors are used
if args.multispec:
    MICASENSE_PATH = args.multispec
else:
    # Default is relative to project location: ../multispec/level0_raw/
    MICASENSE_PATH = Path(proj_file).parents[1] / "multispec/level0_raw"

    if not MICASENSE_PATH.is_dir():
        sys.exit("%s directory does not exist. Check and input paths using -multispec " % str(MICASENSE_PATH))
    else:
        MICASENSE_PATH = str(MICASENSE_PATH)

if args.drtk is not None:
    DRTK_TXT_FILE = args.drtk
    if not Path(DRTK_TXT_FILE).is_file():
        sys.exit("%s file does not exist. Check and input correct path using -drtk option" % str(DRTK_TXT_FILE))

if args.smooth not in DICT_SMOOTH_STRENGTH:
    sys.exit("Value for -smooth must be one of low, medium or high.")

# Export blockshifted P1 positions. Not used in script. Useful for debug or to restart parts of script following any issues.
P1_CAM_CSV = Path(proj_file).parent / "dbg_shifted_p1_pos.csv"
# By default save the CSV with updated MicaSense positions in the MicaSense folder. CSV used within script.
MICASENSE_CAM_CSV = Path(proj_file).parent / "interpolated_micasense_pos.csv"

##################
# Add images
##################
#
# rgb
p1_images = find_files(MRK_PATH, (".jpg", ".jpeg", ".tif", ".tiff"))
chunk = doc.addChunk()
chunk.label = CHUNK_RGB
chunk.addPhotos(p1_images)

# Check that chunk is not empty and images are in default WGS84 CRS
if len(chunk.cameras) == 0:
    sys.exit("Chunk rgb empty")
# check chunk coordinate systems are default EPSG::4326
if "EPSG::4326" not in str(chunk.crs):
    sys.exit("Chunk rgb: script expects images loaded to be in CRS WGS84 EPSG::4326")

#
# multispec
#
micasense_images = find_files(MICASENSE_PATH, (".jpg", ".jpeg", ".tif", ".tiff"))

chunk = doc.addChunk()
chunk.label = CHUNK_MULTISPEC
chunk.addPhotos(micasense_images)
doc.save()

# Check that chunk is not empty and images are in default WGS84 CRS
if len(chunk.cameras) == 0:
    sys.exit("Multispec chunk empty")
if "EPSG::4326" not in str(chunk.crs):
    sys.exit("Multispec chunk: script expects images loaded to be in CRS WGS84 EPSG::4326")

# Check that lever-arm offsets are non-zero:
# As this script is for RGB and MS images captured simultaneously on dual gimbal, lever-arm offsets cannot be 0.
#  Zenmuse P1
if P1_GIMBAL1_OFFSET == 0:
    err_msg = "Lever-arm offset for P1 in dual gimbal mode cannot be 0. Update offset_dict and rerun_script."
    Metashape.app.messageBox(err_msg)

# MicaSense: get Camera Model from one of the images to check the lever-arm offsets for the relevant model
sample_img = open(micasense_images[0], 'rb')
exif_tags = exifread.process_file(sample_img)
cam_model = str(exif_tags.get('Image Model'))

# HARDCODED number of bands.
if len(chunk.sensors) >= 10:
    # Dual sensor (RedEdge-MX Dual: 10, RedEdge-P Dual: 11)
    # Dual sensor: If offsets are 0, exit with error.
    if offset_dict[cam_model]['Dual'] == (0, 0, 0):
        err_msg = "Lever-arm offsets for " + cam_model + " Dual on gimbal 2 cannot be 0. Update offset_dict and rerun script."
        Metashape.app.messageBox(err_msg)
    else:
        MS_GIMBAL2_OFFSET = offset_dict[cam_model]['Dual']
else:
    # RedEdge-MX or RedEdge-P (5-band or 6-band respectively): If offsets are 0, exit with error.
    if offset_dict[cam_model]['Red'] == (0, 0, 0):
        err_msg = "Lever-arm offsets for " + cam_model + " Red on gimbal 2 cannot be 0. Update offset_dict and rerun script."
        Metashape.app.messageBox(err_msg)
    else:
        MS_GIMBAL2_OFFSET = offset_dict[cam_model]['Red']


# Used to find chunks in proc_*
check_chunk_list = [CHUNK_RGB, CHUNK_MULTISPEC]
dict_chunks = {}
for get_chunk in doc.chunks:
    dict_chunks.update({get_chunk.label: get_chunk.key})

# Delete 'Chunk 1' that is created by default.
if 'Chunk 1' in dict_chunks:
    chunk = doc.findChunk(dict_chunks['Chunk 1'])
    doc.remove(chunk)
    doc.save()
#
# # Stop script here. User to click 'Resume Processing' once following steps are complete.
# # see resume_proc() for processing steps.
doc.save()
print("Add images completed.")
print("###########################")
print("###########################")
print("###########################")
print("###########################")
print(
    "Step 1. In the Workspace pane, select multispec chunk. Select Tools-Calibrate Reflectance and 'Locate panels'. Press Cancel once the panels have been located.")
print(
    "Note: The csv of the calibration panel will have to be loaded if this is the first run on the machine. See the protocol for more information.")
print(
    "Step 2. In the Workspace pane under multispec chunk open Calibration images folder. Select and remove images not to be used for calibration.")
print("Step 3. Press the 'Show Masks' icon in the toolbar and inspect the masks on calibration images.")
print(
    "Complete Steps 1 to 3 and press 'Resume Processing' to continue. Reflectance calibration will be completed in the script.")
print("###########################")
print("###########################")
print("###########################")
print("###########################")

label = "Resume processing"
Metashape.app.removeMenuItem(label)
Metashape.app.addMenuItem(label, resume_proc)
Metashape.app.messageBox(
    "Complete Steps 1 to 3 listed on the Console tab and then click on 'Resume Processing' in the toolbar")
