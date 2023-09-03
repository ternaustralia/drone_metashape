# -*- coding: utf-8 -*-
"""
Created May 2023

@author: Poornima Sivanandam

Script to blockshift Metashape chunk coordinates. 
Assumption that images in chunk are in EPSG: 4326 WGS84 coordinate system

Required Input:
    -chunk "<name of Metashape chunk with images to be blockshifted>"
    Example: -chunk "Chunk 1"

    -drtk "<path to file>"
    Path to .txt file containing DRTK field and AUSPOS cartesian coords. 
    Please see Appendix 3 in Drone RGB and Multispectral Imagery Processing Protocol for format of txt file.

Optional Inputs:
    -crs "<EPSG code for target projected coordinate reference system>"
    Example: -crs "7855"


"""

import argparse
import math
import collections
import Metashape
import sys
from pathlib import Path

###############################################################################
# Constants
###############################################################################
GEOG_COORD = collections.namedtuple('GDA2020_GCS', ['lat_decdeg', 'lon_decdeg', 'elliph'])
SOURCE_CRS = Metashape.CoordinateSystem("EPSG::4326")  # WGS84
CONST_a = 6378137  # Semi major axis
#CONST_inv_f = 298.257222101  # Inverse flattening 1/f
CONST_inv_f = 298.257223563  # Inverse flattening 1/f WGS84 ellipsoid



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


def proc_blockshift():
    """
    Author: Poornima Sivanandam
    Arguments: None
    Return: None
    Create: Update coordinates of images in chunk
    Summary:
        * blockshift
        * convert to target CRS (optional through args)
    """
    # Image positions are to be blockshifted, do the following:
    # - Read the .txt file and convert Cartesian coordinates to GDA2020 Lat/Lon
    # - Calculate the difference and apply the shift directly to the cameras (Lon/Lat/Ellipsoidal height) 
    #  (Reference ellipsoids differ for WGS84 and GDA2020 (GRS 1980) - difference in the flattening which results in the semi-minor axis
    #   being different by 0.0001 meters. This is ignored here in applying blockshift calculated using GDA2020 lat/Lon to WGS84 coordinates.)
    # Convert coordinate system for Lat/Lon to target projected coordinate system

    CHUNK_NAME = args.chunk
    chunk = doc.findChunk(dict_chunks[CHUNK_NAME])

    DRTK_TXT_FILE = args.drtk

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
    pos_shift = Metashape.Vector((diff_lon, diff_lat, diff_elliph))

    print("Shifting cameras by: " + str(pos_shift))

    # shift coordinates in the chunk
    for camera in chunk.cameras:
        if not camera.label == camera.master.label:
            continue
        if camera.reference.location:
            camera.reference.location = camera.reference.location + pos_shift

    if args.crs is not None:
        # Convert to projected coordinate system
        target_crs = Metashape.CoordinateSystem("EPSG::" + args.crs)
        for camera in chunk.cameras:
            if not camera.reference.location:
                continue
            camera.reference.location = Metashape.CoordinateSystem.transform(camera.reference.location, SOURCE_CRS,
                                                                             target_crs)
        chunk.crs = target_crs

    # Export updated positions as csv for debug purposes. File not used within script.
    chunk.exportReference(path=str(UPD_POS_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                              delimiter=",", items=Metashape.ReferenceItemsCameras)
    doc.save()

    print("Chunk processing complete!")


############################################
##  Main code
############################################
print("Script start")

# Parse arguments and initialise variables
parser = argparse.ArgumentParser(
    description='Update camera positions in Metashape chunk')
parser.add_argument('-drtk', help='If RGB coordinates to be blockshifted, file containing \
                                                  DRTK base station coordinates from field and AUSPOS', required=True)
parser.add_argument('-chunk', help='Name of chunk to be blockshifted', required=True)
parser.add_argument('-crs',
                    help='EPSG code for target projected CRS. E.g: 7855 for GDA2020/MGA zone 55')

global args
args = parser.parse_args()

# Metashape project
global doc
doc = Metashape.app.document
proj_file = doc.path

DRTK_TXT_FILE = args.drtk
if not Path(DRTK_TXT_FILE).is_file():
    sys.exit("%s file does not exist. Check and input correct path using -drtk option" % str(DRTK_TXT_FILE))

# Export blockshifted positions. Not used in script. Useful for debug or to restart parts of script following any issues.
UPD_POS_CSV = Path(proj_file).parent / "dbg_shifted_pos.csv"

dict_chunks = {}
for get_chunk in doc.chunks:
    dict_chunks.update({get_chunk.label: get_chunk.key})

chunk = doc.findChunk(dict_chunks[args.chunk])
if "EPSG::4326" not in str(chunk.crs):
    sys.exit("Script expects images loaded to be in CRS WGS84 EPSG::4326")

proc_blockshift()

print("End of script")