# -*- coding: utf-8 -*-
"""
Created August 2021

@author: Poornima Sivanandam

Script to process DJI Zenmuse P1 (let gimbal in dual mount gimbal - see GPS/INS offset) .

Assumption that folder structure is as per the TERN protocols:
Data |	Path | Example
Raw data |	<plot>/YYYYMMDD/imagery/<sensor>/level0_raw/ |	SASMDD0001/20220519/imagery/rgb/level0_raw
Data products |	<plot>/YYYYMMDD/imagery/<sensor>/level1_proc/	| SASMDD0001/20220519/imagery/multispec/level1_proc
Metashape project |	plot/YYYYMMDD/imagery/metashape| SASRIV0001/20220516/imagery/metashape/
DRTK logs | plot/YYYYMMDD/drtk/

Raw data paths can be overriden using 'Optional Inputs'.

Required Input:
    -crs "<EPSG code for target projected coordinate reference system>"
    Example: -crs "7855"

Optional Inputs:
    1. -rgb "path to RGB level0_raw folder which also has the MRK file(s)"
        Default is relative to project location: ../rgb/level0_raw/
    2. -smooth "<low/medium/high>"
        Strength value to smooth RGB model. Default is low.
        Low: for low-lying vegetation (grasslands, shrublands), Medium and high: as appropriate for forested sites.
    3. When P1 (RGB camera) coordinates have to be blockshifted:
        - Path to file containing DRTK init and AUSPOS cartesian coords passed using "-drtk <path to file>".

"""

import argparse
import math
import collections
import Metashape
import os
import sys
#from upd_micasense_pos import ret_micasense_pos
from pathlib import Path

# Note: External modules imported were installed through:
# "C:\Program Files\Agisoft\Metashape Pro\python\python.exe" -m pip install <modulename>
# See M300 data processing protocol for more information.

# Metashape Python API updates in v2.0
METASHAPE_V2_PLUS = False
found_version = Metashape.app.version.split('.') # e.g. 2.0.1
if int(found_version[0]) >= 2:
    METASHAPE_V2_PLUS = True

###############################################################################
# Constants
###############################################################################
GEOG_COORD = collections.namedtuple('GDA2020_GCS', ['lat_decdeg', 'lon_decdeg', 'elliph'])

SOURCE_CRS = Metashape.CoordinateSystem("EPSG::4326")  # WGS84

CONST_a = 6378137  # Semi major axis
#CONST_inv_f = 298.257222101  # Inverse flattening 1/f
CONST_inv_f = 298.257223563  # Inverse flattening 1/f WGS84 ellipsoid

# Chunks in Metashape
CHUNK_RGB = "rgb"

IMG_QUAL_THRESHOLD = 0.7

DICT_SMOOTH_STRENGTH = {'low': 50, 'medium': 100, 'high': 200}


###############################################################################
# Function definitions
###############################################################################
def cartesian_to_geog(X, Y, Z):
    """
    Author: Poornima Sivanandam
    Convert Cartesian coordinates to GDA2020 and return Lat, Lon, ellipsoidal height as a named tuple.
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
    # - Read the .txt file and convert Cartesian coordinates to GDA2020 Lat/Lon
    # - Calculate the difference and apply the shift directly to the cameras (Lon/Lat/Ellipsoidal height) in 'P1' chunk
    #  (Reference ellipsoids differ for WGS84 and GDA2020 (GRS 1980) - difference in the flattening which results in the semi-minor axis
    #   being different by 0.0001 meters. This is ignored here in applying blockshift calculated using GDA2020 lat/Lon to WGS84 coordinates.)
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
        diff_lat = round((drtk_auspos.lat_decdeg - drtk_field.lat_decdeg), 8)
        diff_lon = round((drtk_auspos.lon_decdeg - drtk_field.lon_decdeg), 8)
        diff_elliph = round((drtk_auspos.elliph - drtk_field.elliph), 6)
        P1_shift = Metashape.Vector((diff_lon, diff_lat, diff_elliph))

        print("Shifting P1 cameras by: " + str(P1_shift))

        # shift coordinates in the chunk
        for camera in chunk.cameras:
            if not camera.label == camera.master.label:
                continue
            if camera.reference.location:
                camera.reference.location = camera.reference.location + P1_shift

    # Convert to projected coordinate system
    target_crs = Metashape.CoordinateSystem("EPSG::" + args.crs)
    for camera in chunk.cameras:
        if not camera.reference.location:
            continue
        camera.reference.location = Metashape.CoordinateSystem.transform(camera.reference.location, SOURCE_CRS,
                                                                         target_crs)

    chunk.crs = target_crs

    if blockshift_p1:
        # Export updated positions as csv for debug purposes. Not used in script.
        chunk.exportReference(path=str(P1_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                              delimiter=",", items=Metashape.ReferenceItemsCameras)
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
    chunk.sensors[0].antenna.location_ref = Metashape.Vector((0.087, 0.0, 0.0))
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
    # check if exist and reuse depthmap? # reuse_depth=True below
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
    chunk.buildModel(surface_type=Metashape.HeightField, face_count=Metashape.MediumFaceCount)
    doc.save()

    # Decimate and smooth mesh to use as orthorectification surface
    # Halve face count?
    chunk.decimateModel(face_count=len(chunk.model.faces) / 2)
    # Smooth model
    smooth_val = DICT_SMOOTH_STRENGTH[args.smooth]
    chunk.smoothModel(smooth_val)

    #
    # Build and export orthomosaic
    #
    print("Build orthomosaic")
    chunk.buildOrthomosaic(surface_data=Metashape.DataSource.ModelData, refine_seamlines=True)
    doc.save()

    if chunk.orthomosaic:
        # Round resolution to 2 decimal places
        res_xy = round(chunk.orthomosaic.resolution, 2)

        # if p1/ folder does not exist in MRK_PATH save orthomosaic in the project directory
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


############################################
##  Main code
############################################
print("Script start")

# Parse arguments and initialise variables
parser = argparse.ArgumentParser(
    description='Update camera positions in P1 rgb chunk in Metashape project')
parser.add_argument('-crs',
                    help='EPSG code for target projected CRS. E.g: 7855 for GDA2020/MGA zone 55',
                    required=True)
parser.add_argument('-rgb', help='path to RGB level0_raw folder that also has the MRK files')
parser.add_argument('-smooth', help='Smoothing strength used to smooth RGB mesh low/med/high', default="low")
parser.add_argument('-drtk', help='If RGB coordinates to be blockshifted, file containing \
                                                  DRTK base station coordinates from field and AUSPOS')

global args
args = parser.parse_args()
global MRK_PATH

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

# HARDCODED sensor names
if args.rgb:
    MRK_PATH = args.rgb
else:
    # Default is relative to project location: ../rgb/level0_raw/
    MRK_PATH = Path(proj_file).parents[1] / "rgb/level0_raw"
    if not MRK_PATH.is_dir():
        sys.exit("%s directory does not exist. Check and input paths using -rgb " % str(MRK_PATH))
    else:
        MRK_PATH = str(MRK_PATH)

if args.drtk is not None:
    DRTK_TXT_FILE = args.drtk
    if not Path(DRTK_TXT_FILE).is_file():
        sys.exit("%s file does not exist. Check and input correct path using -drtk option" % str(DRTK_TXT_FILE))

if args.smooth not in DICT_SMOOTH_STRENGTH:
    sys.exit("Value for -smooth must be one of low, medium or high.")

# Export blockshifted P1 positions. Not used in script. Useful for debug or to restart parts of script following any issues.
P1_CAM_CSV = Path(proj_file).parent / "dbg_shifted_p1_pos.csv"


##################
# Add images
##################
#
# p1
#
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

doc.save()

# Used to find chunks in proc_*
check_chunk_list = [CHUNK_RGB]
dict_chunks = {}
for get_chunk in doc.chunks:
    dict_chunks.update({get_chunk.label: get_chunk.key})

# Delete 'Chunk 1' that is created by default.
if 'Chunk 1' in dict_chunks:
    chunk = doc.findChunk(dict_chunks['Chunk 1'])
    doc.remove(chunk)
    doc.save()

print("Add images completed.")

proc_rgb()

print("End of script")