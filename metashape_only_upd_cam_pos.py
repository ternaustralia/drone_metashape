# -*- coding: utf-8 -*-
"""
Created August 2021

@author: Poornima Sivanandam

Script to update camera positions in chunks p1 and micasense in the Metashape project.

Assumption that folder structure is as per the TERN protocols:
Data |	Path | Example
Raw data |	<plot>/YYYYMMDD/imagery/<sensor>/level0_raw/ |	SASMDD0001/20220519/imagery/rgb/level0_raw
Data products |	<plot>/YYYYMMDD/imagery/<sensor>/level1_proc/	| SASMDD0001/20220519/imagery/multispec/level1_proc
Metashape project |	plot/YYYYMMDD/imagery/metashape| SASRIV0001/20220516/imagery/metashape/

Raw data paths can be overriden using 'Optional Inputs'.

Required Input:
    -crs "<EPSG code for target projected coordinate reference system - used in MicaSense position interpolation>"
    Example: -crs "7855"

Optional Inputs:
    1. -multispec "path to multispectral level0_raw folder containing raw data"
        Default is relative to project location: ../multispec/level0_raw/
    2. -rgb "path to RGB level0_raw folder which also has the MRK file(s)"
        Default is relative to project location: ../rgb/level0_raw/
    3. When P1 coordinates have to be blockshifted:
        - path to file containing DRTK init and AUSPOS cartesian coords passed using "-drtk <path to txt file>".

Summary:
Note that following steps must be done before running script
    1. Add RGB and multispectral images.
    2. Micasense - Locate panels, Calibration images selected and Reflectance calibration complete.

chunk rgb:
When '-drtk <path to txt file>' is passed:
    - In the .txt file, the first line must have the comma separated cartesian coordinates from RINEX OBS file and the
    second line the coordinates from AUSPOS report. Contents of .txt file as an example :
        -3943357.3162, 2612482.7744, -4265086.5255
        -3943357.844, 2612483.133, -4265086.641
    - The assumption is that the image coordinates are in Lat/Lon - the default geographic coordinate system.
    - P1 Camera positions are then updated directly in the Metashape project. 
CRS updated to the target Coordinate Reference System set through target_crs
   
chunk multispec:
    - MicaSense imagery are geotagged with navigation-grade accuracy. This script will interpolate positions of the MicaSense sensor using the P1 MRK files.
    - We use the MRK files rather than EXIF tags from the P1 images (or the blockshifted images in the Metashape project) because the P1 EXIF tags 
    do not contain sub-second time that is required for interpolation.
    - However as MRK files contain positions recorded in the field, the blockshift (if enabled) has to be added separately to these 
    coordinates. This is done through an argument to ret_micasense_pos.
    - The master band in Micasense chunk is identified. Image files corresponding to the master band suffix are read from the 'micasense'
    path which is the first positional argument to this script.
    In module inteprolate_micasense_coord (function ret_micasense_pos):
        - EXIF tags of these image files are read for TimeStamp information
        - P1 posiitons and timestamps are read from P1 MRK
        - Camera positions are converted to the target projected coordinate system - necessary for interpolation based on timestamp
        - Micasense position calculated using the two closest (beofre and after) P1 images.
        - Blockshift (if enabled) applied to Micasense position.
        - MicaSense capture is triggered separately, hence there will be images that triggered when P1 did not. Interpolation is not possible 
        for these images. As our goal is co-registered P1 and MicaSense imagery, these images can be deleted. This is done as follows:
            - if MicaSense image timestamps are outside P1 times, a dummy 'Altitude' value of 0 is added. This value is used in the main script to
            delete these images
            - This is also why images used for refelctance calibration are required to first be identified and moved to the 'Calibration images' 
            group before this script is run (see M300 Data processing protocol). Otherwise as these images triggered outside P1 times, they will be deleted.
    - Finally ret_micasense_pos returns a CSV file with the updated Easting, Northing, Altitude for all MicaSense master band images.
    - The updated positions are loaded in the micasense chunk
    - Images with Altitude 0 triggered outside P1 times and are hence deleted.


"""

import argparse
import math
import collections
import numpy as np
import Metashape
import sys
from pathlib import Path
from upd_micasense_pos import ret_micasense_pos

# Note: External modules imported were installed through:
# "C:\Program Files\Agisoft\Metashape Pro\python\python.exe" -m pip install <modulename> 
# See M300_data_processing protocol for more information.

###############################################################################
# Constants
###############################################################################
GEOG_COORD = collections.namedtuple('GDA2020_GCS', ['lat_decdeg', 'lon_decdeg', 'elliph'])

SOURCE_CRS = Metashape.CoordinateSystem("EPSG::4326")  # WGS84
GDA2020_CRS = Metashape.CoordinateSystem("EPSG::7844")  # GDA2020

CONST_a = 6378137  # Semi major axis
#CONST_inv_f = 298.257222101  # Inverse flattening 1/f
CONST_inv_f = 298.257223563  # Inverse flattening 1/f WGS84 ellipsoid

# Edit these values if necessary to refer to the correct P1/Micasense chunks in the project.
CHUNK_MULTISPEC = "multispec"
CHUNK_RGB = "rgb"


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


############################################
##  Main code
############################################
print("Script start")
parser = argparse.ArgumentParser(
    description='Update camera positions in P1 and/or MicaSense chunks in Metashape project')
parser.add_argument('-crs',
                    help='EPSG code for target projected CRS for micasense cameras. E.g: 7855 for GDA2020/MGA zone 55',
                    required=True)
parser.add_argument('-drtk', help='If RGB coordinates to be blockshifted, file containing \
                                                  DRTK base station coordinates from field and AUSPOS')

args = parser.parse_args()
TARGET_CRS = Metashape.CoordinateSystem("EPSG::" + args.crs)

if args.drtk is not None:
    blockshift_p1 = True
    DRTK_TXT_FILE = args.drtk
    if not Path(DRTK_TXT_FILE).is_file():
        sys.exit("%s file does not exist. Check and input correct path using -drtk option" % str(DRTK_TXT_FILE))
else:
    blockshift_p1 = False

# Metashape project
doc = Metashape.app.document
proj_file = doc.path

# Export blockshifted P1 positions. Not used in script. Useful for debug or to restart parts of script following any issues.
P1_CAM_CSV = Path(proj_file).parent / "shifted_p1_pos.csv"

# By default save the CSV with updated MicaSense positions in the MicaSense folder. CSV used within script.
MICASENSE_CAM_CSV = Path(proj_file).parent / "interpolated_micasense_pos.csv"

#####################
# Checks
#####################
check_chunk_list = [CHUNK_RGB, CHUNK_MULTISPEC]
dict_chunks = {}

for get_chunk in doc.chunks:
    dict_chunks.update({get_chunk.label: get_chunk.key})

# Delete 'Chunk 1' that is created by default.
if 'Chunk 1' in dict_chunks:
    chunk = doc.findChunk(dict_chunks['Chunk 1'])
    doc.remove(chunk)
    doc.save()

# If rgb and multispec chunks not present, exit script.
chunks_exist = all(item in dict_chunks for item in check_chunk_list)
if not chunks_exist:
    sys.exit(
        "Ensure that rgb and multispectral images have been added and the chunks renamed as rgb and multispec respectively.")

# Check that images have been added
chunk = doc.findChunk(dict_chunks[CHUNK_RGB])
if len(chunk.cameras) == 0:
    sys.exit("Chunk rgb empty")
# check chunk coordinate systems are default EPSG::4326
if "EPSG::4326" not in str(chunk.crs):
    sys.exit("Chunk p1: script expects images loaded to be in CRS WGS84 EPSG::4326")

# Repeat for micasense
chunk = doc.findChunk(dict_chunks[CHUNK_MULTISPEC])
if len(chunk.cameras) == 0:
    sys.exit("Chunk multispec empty")
if "EPSG::4326" not in str(chunk.crs):
    sys.exit("Chunk multispec: script expects images loaded to be in CRS WGS84 EPSG::4326")

######################
##  P1
######################
# If P1 positions are to be blockshifted, do the following:
# - Duplicate chunk 'P1'
# - Read the .txt file and convert Cartesian coordinates to GDA2020 Lat/Lon
# - Calculate the difference and apply the shift directly to the cameras (Lon/Lat/Ellipsoidal height) in 'P1' chunk
#  (Reference ellipsoids differ for WGS84 and GDA2020 (GRS 1980) - difference in the flattening which results in the semi-minor axis
#   being different by 0.0001 meters. This is ignored here in applying blockshift calculated using GDA2020 lat/Lon to WGS84 coordinates.)
# Convert coordinate system for Lat/Lon to target projected coordinate system (not necessary - done for consistency with micasense chunk).
#
chunk = doc.findChunk(dict_chunks[CHUNK_RGB])

# Get path of RGB images
rgb_image = chunk.cameras[0]
rgb_image_path = Path(rgb_image.photo.path)
MRK_PATH = str(rgb_image_path.parents[1]) #  HARDCODED path to RGB level_0 folder

if blockshift_p1:
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
for camera in chunk.cameras:
    if not camera.reference.location:
        continue
    camera.reference.location = Metashape.CoordinateSystem.transform(camera.reference.location, SOURCE_CRS, TARGET_CRS)

chunk.crs = TARGET_CRS

if blockshift_p1:
    # Export updated positions as csv for debug purposes. Not used in script.
    chunk.exportReference(path=str(P1_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                          delimiter=",", items=Metashape.ReferenceItemsCameras)
doc.save()

######################
##  MicaSense
######################
# Move to micasense chunk
chunk = doc.findChunk(dict_chunks[CHUNK_MULTISPEC])

# Get path of multispec images
multispec_image = chunk.cameras[0]
multispec_image_path = Path(multispec_image.photo.path)
MICASENSE_PATH = str(multispec_image_path.parents[2]) # HARDCODED path to the multispec level_0 folder

# Get image suffix of master camera
camera = chunk.cameras[0]
cam_master = camera.master.label.split('_')

# Image file naming assumption: IMG_xxxx_suffixNum
img_suffix_master = cam_master[2]
print(img_suffix_master)

# If P1 coords were block shifted, pass vector for x, y, z shift
# else vector is 0
if blockshift_p1:
    P1_shift_vec = np.array([diff_lat, diff_lon, diff_elliph])
else:
    P1_shift_vec = np.array([0.0, 0.0, 0.0])

print("Interpolate Micasense position based on P1")
print("With blockshift for P1: " + str(P1_shift_vec))

# inputs: paths to MRK for P1 position, Micasense images on the drive, image suffix to get master band image files, target CRS
# returns output csv file with updated micasense positions
ret_micasense_pos(MRK_PATH, MICASENSE_PATH, img_suffix_master, args.crs,
                  str(MICASENSE_CAM_CSV), P1_shift_vec)

# Load updated positions in the chunk
chunk.importReference(str(MICASENSE_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                      delimiter=",", crs=TARGET_CRS, skip_rows=1,
                      items=Metashape.ReferenceItemsCameras)
doc.save()

# ret_micasense_pos wrote Altitude = 0 (last column) for MicaSense cams that triggered when P1 did not.
# Get list of cameras with Altitude = 0
del_camera_names = list()

# Only look at altitude of master band images
for camera in chunk.cameras:
    if not camera.label == camera.master.label:
        continue
    if not camera.reference.location:
        continue
    if camera.reference.location.z == 0:
        del_camera_names.append(camera.label)

# Delete the images but do not remove any calibration images
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

print("End of script. Updated micasense camera positions have been imported in Reference pane.")
