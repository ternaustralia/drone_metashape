# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 14:15:33 2021

@author:  Poornima Sivanandam


NOTE:   This is a prototype on creating custom widgets within Metashape.
        Comments have not been added yet.
        Code needs re-organising - currently copied functions and code as is from the other metashape_* scripts.

To run the script, copy this script and up_micasense_pos.py to AppData/Local/Agisoft/Metashape Pro/scripts/ folder under
your username folder on the C drive.

Follow steps per M300_data_processing document.

"""

import Metashape
from PySide2 import QtCore, QtWidgets
from PySide2.QtWidgets import *
import os
import glob
import sys
from pathlib import Path
import math
import collections
import numpy as np
from upd_micasense_pos import ret_micasense_pos


# # TODO auto check compatibility based on latest version tested
# compatible_major_version = "1.7"
# found_major_version = ".".join(Metashape.app.version.split('.')[:2])
# if found_major_version != compatible_major_version:
#     raise Exception("Incompatible Metashape version: {} != {}".format(found_major_version, compatible_major_version))

global path_tern_plot
###############################################################################
# Constants
###############################################################################
GDA2020_COORD = collections.namedtuple('GDA2020_GCS', ['lat_decdeg', 'lon_decdeg', 'elliph'])

SOURCE_CRS = Metashape.CoordinateSystem("EPSG::4326")  # WGS84

CONST_a = 6378137  # Semi major axis
CONST_inv_f = 298.257222101  # Inverse flattening 1/f

IMG_QUAL_THRESHOLD = 0.7

# Edit these values if necessary to refer to the correct P1/Micasense chunks in the project.
CHUNK_RGB = "rgb"
CHUNK_MULTISPEC = "multispec"

DICT_SMOOTH_STRENGTH = {'Low': 50, 'Medium': 100, 'High': 200}


###############################################################################
# Function definitions
###############################################################################
def cartesian_to_gda2020(X, Y, Z):
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

    conv_coord = GDA2020_COORD(lat_dec_deg, lon_dec_deg, ellip_h)

    return conv_coord


def startProc():

    global doc
    doc = Metashape.app.document
    app = QApplication.instance()
    window = mainMetashapeWindow()
    window.show()


class addChunksDlg(QDialog):
    def __init__(self):
        super().__init__()

        labelMenu = QLabel()
        labelMenu.setText("Adding images")
        clickAddImages = QPushButton("Create chunks")
        clickQuit = QPushButton("Cancel")

        layout = QGridLayout()
        widgetList = [labelMenu, clickAddImages, clickQuit]
        for w in widgetList:
            layout.addWidget(w)
        self.setLayout(layout)

        clickAddImages.clicked.connect(lambda: self.processaddImages())
        clickQuit.clicked.connect(lambda: self.close())
        self.exec()

    def processaddImages(self):
        # add chunks
        global path_tern_plot
        path_tern_plot = QtWidgets.QFileDialog.getExistingDirectory(self,
                                                                    'Choose the folder where rgb/ and multispec/ folders are located.')

        p1_path = os.path.join(path_tern_plot, CHUNK_RGB)
        if not os.path.isdir(p1_path):
            sys.exit("Cannot find rgb/ in TERN plot folder. Check folder location and naming.")
        p1_images = self.find_files(p1_path, (".jpg", ".jpeg", ".tif", ".tiff"))
        chunk = doc.addChunk()
        chunk.label = CHUNK_RGB
        chunk.addPhotos(p1_images)

        micasense_path = os.path.join(path_tern_plot, CHUNK_MULTISPEC)
        if not os.path.isdir(micasense_path):
                sys.exit(
                    "Cannot find multispec/ in TERN plot folder. Check folder location and naming.")

        micasense_images = self.find_files(micasense_path, (".jpg", ".jpeg", ".tif", ".tiff"))
        chunk = doc.addChunk()
        chunk.label = CHUNK_MULTISPEC
        chunk.addPhotos(micasense_images)

        proj_path = os.path.join(path_tern_plot, "metashape/")
        if not os.path.exists(proj_path):
            os.makedirs(proj_path)
        proj_name = Metashape.app.getSaveFileName("Specify the project name:", dir=proj_path)
        doc.save(proj_name)
        print("Add images completed.")
        print("Complete the other steps as listed on the 'Metashape processing menu'.")
        self.close()

    def find_files(self, folder, types):
        photo_list = list()
        for dir, subdir, file in os.walk(folder):
            for filename in file:
                if (filename.lower().endswith(types)):
                    photo_list.append(os.path.join(dir, filename))
        return (photo_list)


class promptReflCalibDlg(QWidget):
    def __init__(self):
        super().__init__()

        label = QLabel()
        label.setText("Steps to complete using Metashape toolbar")

        msg1 = QLabel(
            "Step 1: In the Workspace pane, select multispec chunk. Select Tools-Calibrate Relectance and 'Locate panels'. Once this is done, press Cancel.")
        msg2 = QLabel(
            "Step 2: In the Workspace pane under multispec chunk open Calibration images folder. Select and remove images not to be used for calibration.")
        msg3 = QLabel(
            "Step 3: Press the 'Show Masks' icon in the tool bar and inspect the masks on calibration images.")
        msg4 = QLabel(
            "Reflectance calibration will be completed in the script. Enter other inputs in the Metashape processing menu and click 'Run'")

        clickOK = QPushButton("OK")

        layout = QGridLayout()
        widgetList = [label, msg1, msg2, msg3, msg4, clickOK]

        for w in widgetList:
            layout.addWidget(w)

        self.setLayout(layout)
        clickOK.clicked.connect(lambda: self.close())


class mainMetashapeWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Metashape processing menu")

        layout = QVBoxLayout()

        clickAddImages = QPushButton("Add images")
        self.label = QLabel("User input required to Locate panel and select calibration images")
        info = QPushButton("Info on steps to complete")
        self.chkLocatePanel = QCheckBox("Tools->Calibrate Reflectance - Locate panels complete?")
        self.chkCalibImages = QCheckBox("Calibration Images inspected and images removed as needed?")

        # # type of MicaSense sensor
        # self.labelSensor = QLabel()
        # self.labelSensor.setText("Select Micasense Sensor")
        # self.micasenseSensor = QComboBox()
        # micasenseSensorOpts = ["Dual", "RedEdgeMX"]
        # for val in micasenseSensorOpts:
        #     self.micasenseSensor.addItem(val)
        # self.micasenseSensor.setCurrentText("RedEdgeMX")

        # target CRS
        self.labelCRS = QLabel()
        self.labelCRS.setText("Enter EPSG code of output projected coordinate system")
        self.targetCRS = QLineEdit()
        self.targetCRS.setPlaceholderText(
            "EPSG code of output projected coordinate system. E.g.7855 for GDA2020 MGAZone 55")
        self.targetCRS.setMaxLength(4)

        # smooth strength
        self.labelSmooth = QLabel()
        self.labelSmooth.setText("Select Strength option used to smooth P1 model")
        self.smoothStrength = QComboBox()
        smoothStrengthOpts = ["Low", "Medium", "High"]
        for val in smoothStrengthOpts:
            self.smoothStrength.addItem(val)
        self.smoothStrength.setCurrentText("Low")

        self.blockshiftP1 = QCheckBox("Blockshift P1 image coordinates using AUSPOS results?")
        self.labelDRTKFile = QLabel(
            "If blockshift enabled, choose txt file containing DRTK Cartesian coords from field and AUSPOS")
        self.labelCRS.setText("Enter EPSG code of output projected coordinate system")
        # drtkTxtFile = QtWidgets.QFileDialog.getOpenFileFileName(self, '(If Blockshift enabled) Choose text file containing DRTK cartesian coords from field and AUSPOS',
        #                                                              filter="Text file (*.txt)")
        self.drtkTxtFile = QLineEdit()
        self.drtkTxtFile.setPlaceholderText('')
        openDrtkTxtFile = QPushButton("Browse")


        clickStart = QPushButton("Run")
        clickQuit = QPushButton("Quit")

        widgetList = [clickAddImages, self.label, info, self.chkLocatePanel, self.chkCalibImages,
                      self.labelCRS, self.targetCRS,
                      self.labelSmooth, self.smoothStrength, self.blockshiftP1, self.labelDRTKFile, self.drtkTxtFile,
                      openDrtkTxtFile, clickStart, clickQuit]

        for w in widgetList:
            layout.addWidget(w)

        widget = QWidget()
        widget.setLayout(layout)

        clickAddImages.clicked.connect(lambda: self.proc_openNewWindow())
        info.clicked.connect(lambda: self.proc_displayInfo())
        openDrtkTxtFile.clicked.connect(lambda: self.openDrtkFileDlg())
        clickStart.clicked.connect(lambda: self.processChunkWorkflow())
        clickQuit.clicked.connect(lambda: self.close())

        self.setCentralWidget(widget)


    def proc_openNewWindow(self):
        self.dlg = addChunksDlg()
        #TODO: call promptReflCalibDlg() here instead of separate Info button?


    def proc_displayInfo(self):
        self.dlg2 = promptReflCalibDlg()
        self.dlg2.show()


    def processChunkWorkflow(self):
        print("chunk workflow")
        self.locatePanelDone = self.chkLocatePanel.isChecked()
        self.calibImageSelDone = self.chkCalibImages.isChecked()

        if self.locatePanelDone & self.calibImageSelDone:
            print("User input for reflectance calibration complete. Resuming image processing.")
        else:
            sys.exit("Locate panel, and select appropriate calibration images.")

        global dict_chunks
        dict_chunks = {}

        for get_chunk in doc.chunks:
            dict_chunks.update({get_chunk.label: get_chunk.key})

        #TODO self.check_chunks()
        self.proc_chunks()


    def check_chunks(self):
        # If rgb and multispec chunks not present, exit script.

        chunks_exist = all(item in dict_chunks for item in [CHUNK_RGB, CHUNK_MULTISPEC])
        if not chunks_exist:
            sys.exit(
                "Ensure that rgb and multispec images have been added and the chunks renamed as rgb and multispec respectively.")

        # Check that chunks are not empty
        chunk = doc.findChunk(dict_chunks[CHUNK_RGB])
        if len(chunk.cameras) == 0:
            sys.exit("Chunk rgb is empty")

        # check chunk coordinate systems are default EPSG::4326
        if "EPSG::4326" not in str(chunk.crs):
            sys.exit("Chunk rgb: script expects images loaded to be in WGS84 EPSG::4326")

        # Repeat for micasense
        chunk = doc.findChunk(dict_chunks[CHUNK_MULTISPEC])
        if len(chunk.cameras) == 0:
            sys.exit("Chunk multispec is empty")

        if "EPSG::4326" not in str(chunk.crs):
            sys.exit("Chunk multispec: script expects images loaded to be in WGS84 EPSG::4326")

        print("check_chunks complete!")


    def openDrtkFileDlg(self):
        filename = QFileDialog.getOpenFileName(self,
                                               '(If Blockshift enabled) Choose text file containing DRTK cartesian coords from field and AUSPOS',
                                               filter="Text file (*.txt)")
        if len(filename[0]) > 3:
            self.drtkTxtFile.setText(filename[0])

    def proc_chunks(self):
        ######################
        ##  P1
        ######################
        chunk = doc.findChunk(dict_chunks[CHUNK_RGB])
        proj_file = doc.path

        blockshift_p1 = self.blockshiftP1.isChecked()
        CARTESIAN_COORDS = self.drtkTxtFile.text()
        EPSG_CRS = self.targetCRS.text()
        TARGET_CRS = Metashape.CoordinateSystem("EPSG::" + EPSG_CRS)

        if blockshift_p1:
            # read from txt/csv cartesian for RTK initial (line 1) and AUSPOS coords (line 2)
            with open(CARTESIAN_COORDS, 'r') as file:
                line = file.readline()
                split_line = line.split(',')
                init_gda2020 = cartesian_to_gda2020(float(split_line[0]), float(split_line[1]), float(split_line[2]))
                line = file.readline()
                split_line = line.split(',')
                auspos_gda2020 = cartesian_to_gda2020(float(split_line[0]), float(split_line[1]), float(split_line[2]))

            # calc difference
            diff_lat = round((auspos_gda2020.lat_decdeg - init_gda2020.lat_decdeg), 6)
            diff_lon = round((auspos_gda2020.lon_decdeg - init_gda2020.lon_decdeg), 6)
            diff_elliph = round((auspos_gda2020.elliph - init_gda2020.elliph), 6)
            P1_shift = Metashape.Vector((diff_lon, diff_lat, diff_elliph))
            # P1_shift = Metashape.Vector([-3e-06, -4e-06, -1.677178])

            print("Shifting P1 cameras by: " + str(P1_shift))

            # shift coordinates in the chunk
            for camera in chunk.cameras:
                if not camera.label == camera.master.label:
                    continue
                if not camera.reference.location:
                    continue
                else:
                    camera.reference.location = camera.reference.location + P1_shift

            # Export updated positions as csv for debug purposes. Not used in script.
            # Export blockshifted P1 positions. Not used in script. Useful for debug or to restart parts of script following any issues.
            P1_CAM_CSV = Path(proj_file).parent / "shifted_p1_pos.csv"

            chunk.exportReference(path=str(P1_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                                  delimiter=",", items=Metashape.ReferenceItemsCameras)

        # Convert to projected coordinate system
        for camera in chunk.cameras:
            if not camera.reference.location:
                continue
            camera.reference.location = Metashape.CoordinateSystem.transform(camera.reference.location, SOURCE_CRS,
                                                                             TARGET_CRS)

        chunk.crs = TARGET_CRS
        doc.save()

        #
        # Estimate image quality and remove cameras with quality < threshold
        #
        chunk.analyzePhotos()
        low_img_qual = [camera for camera in chunk.cameras if
                        (float(camera.meta["Image/Quality"]) < IMG_QUAL_THRESHOLD)]
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
        chunk.buildDenseCloud()
        doc.save()

        #
        # Build Mesh
        #
        print("Build mesh")
        chunk.buildModel(surface_type=Metashape.HeightField, face_count=Metashape.MediumFaceCount)
        doc.save()

        # Decimate and smooth mesh to use as orthorectification surface
        # TODO: is this final? decimate model by half
        chunk.decimateModel(face_count=len(chunk.model.faces) / 2)

        smooth_val = DICT_SMOOTH_STRENGTH[self.smoothStrength.currentText()]
        chunk.smoothModel(smooth_val)

        # Export model for use in micasense chunk
        model_file = Path(proj_file).parent / (Path(proj_file).stem + "_rgb_smooth_" + str(smooth_val) + ".obj")
        chunk.exportModel(path=str(model_file), crs=TARGET_CRS, format=Metashape.ModelFormatOBJ)

        #
        # Build and export orthomosaic
        #
        print("Build orthomosaic")
        #chunk.buildOrthomosaic(surface_data=Metashape.DataSource.ModelData, refine_seamlines=True)
        chunk.buildOrthomosaic(refine_seamlines=True)
        doc.save()

        if chunk.orthomosaic:
            # Round resolution to 2 decimal places
            res_xy = round(chunk.orthomosaic.resolution, 2)

            # create level1_proc folder if it does not exist
            Path(Path(proj_file).parents[1] / "rgb/level1_proc/").mkdir(parents=True, exist_ok=True)
            # file naming format: <projname>_rgb_ortho_<res_in_m>.tif
            ortho_file = Path(proj_file).parents[1] / "rgb/level1_proc/" / (
                    Path(proj_file).stem + "_rgb_ortho_" + str(res_xy).split('.')[1] + ".tif")

            compression = Metashape.ImageCompression()
            compression.tiff_compression = Metashape.ImageCompression.TiffCompressionDeflate
            compression.tiff_big = True
            compression.tiff_tiled = True
            compression.tiff_overviews = True

            chunk.exportRaster(path=str(ortho_file), resolution_x=res_xy, resolution_y=res_xy,
                               image_format=Metashape.ImageFormatTIFF,
                               save_alpha=False, source_data=Metashape.OrthomosaicData, image_compression=compression)
            print("Exported orthomosaic: " + str(ortho_file))
            print("RGB chunk processing complete!")

        ######################
        ##  MicaSense
        ######################

        chunk = doc.findChunk(dict_chunks[CHUNK_MULTISPEC])

        # Get image suffix of master camera
        camera = chunk.cameras[0]
        cam_master = camera.master.label.split('_')

        # file naming assumption: IMG_xxxx_suffixNum
        img_suffix_master = cam_master[2]
        print(img_suffix_master)

        # If P1  blockshifted, pass vector for x, y, z shift
        if blockshift_p1:
            P1_shift_vec = np.array([diff_lat, diff_lon, diff_elliph])
        else:
            P1_shift_vec = np.array([0.0, 0.0, 0.0])

        print("Interpolate Micasense position based on P1 with blockshift" + str(P1_shift_vec))

        # inputs: paths to MRK for P1 position, Micasense images on the drive, image suffix to get master band image files, target CRS
        # returns output csv file with updated micasense positions
        MICASENSE_CAM_CSV = Path(proj_file).parent / "interpolated_micasense_pos.csv"

        # path_tern_plot = Path(proj_file).parents[1]

        MRK_PATH = os.path.join(path_tern_plot, CHUNK_RGB)
        MICASENSE_PATH = os.path.join(path_tern_plot, CHUNK_MULTISPEC)
        ret_micasense_pos(MRK_PATH, MICASENSE_PATH, img_suffix_master, EPSG_CRS,
                          str(MICASENSE_CAM_CSV), P1_shift_vec)

        # Load updated positions in the chunk
        chunk.importReference(str(MICASENSE_CAM_CSV), format=Metashape.ReferenceFormatCSV, columns="nxyz",
                              delimiter=",", crs=TARGET_CRS, skip_rows=1,
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

        #
        # Set primary channel
        #
        # Get index of NIR band. Micasense Dual: NIR is sensors[9], and in RedEdge-M sensors[4]
        for s in chunk.sensors:
            if s.label.find("NIR") != -1:
                nir_idx = s.layer_index
                break

        # channels are numbered from 1 to ... So add 1 to idx for channel.
        print("Setting primary channel to " + s.label)
        chunk.primary_channel = nir_idx + 1

        #
        # GPS/INS offset for master sensor
        #
        # Offsets for RedEdge-MX Dual: (-0.097, 0.02, -0.08)
        #            RedEdge-MX: (-0.097, -0.03, -0.06)
        print("Updating Micasense GPS offset")
        if len(chunk.sensors) == 10:
            chunk.sensors[0].antenna.location_ref = Metashape.Vector((-0.097, 0.02, -0.08))
            # NOT USED - was used in naming ortho tif. Updated to generic "multispec"
            micasense_sensor = "micasense_dual"
        else:
            chunk.sensors[0].antenna.location_ref = Metashape.Vector((-0.097, -0.03, -0.06))
            micasense_sensor = "micasense"

        #
        # Set Raster Transform to calculate reflectance
        #
        print("Updating Raster Transform for relative reflectance")
        raster_transform_formula = []
        num_bands = len(chunk.sensors)
        for band in range(1, num_bands + 1):
            raster_transform_formula.append("B" + str(band) + "/32768")
        chunk.raster_transform.formula = raster_transform_formula
        chunk.raster_transform.calibrateRange()
        chunk.raster_transform.enabled = True
        doc.save()

        #
        # Estimate image quality and remove cameras with quality < threshold
        #

        chunk.analyzePhotos()

        low_img_qual = [camera.master for camera in chunk.cameras if (float(camera.meta["Image/Quality"]) < 0.5)]
        if low_img_qual:
            print("Removing cameras with Image Quality < %.1f" % 0.5)
            chunk.remove(list(set(low_img_qual)))
        doc.save()

        # LocatePanel and inspecting calibration images/masks completed by user prior to this
        # See locatePanelDone and calibImageSelDone
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
        # Build Dense Cloud
        #
        # Downscale: ultra, high, medium, low, lowest: 1, 2, 4, 8, 16
        print("Build dense cloud")
        chunk.buildDepthMaps(downscale=4)  # medium quality. and default: mild filtering.
        chunk.buildDenseCloud()
        doc.save()

        #
        # Build and export orthomosaic
        #
        # Import P1 model for use in orthorectification
        chunk.importModel(path=str(model_file), crs=TARGET_CRS, format=Metashape.ModelFormatOBJ)

        print("Build orthomosaic")
        chunk.buildOrthomosaic(surface_data=Metashape.DataSource.ModelData, refine_seamlines=True)
        doc.save()

        if chunk.orthomosaic:
            # Round resolution to 2 decimal places
            res_xy = round(chunk.orthomosaic.resolution, 2)

            # create level1_proc folder if it does not exist
            Path(Path(proj_file).parents[1] / "multispec/level1_proc/").mkdir(parents=True, exist_ok=True)

            # file naming format: <projname>_multispec_ortho_<res_in_m>.tif
            ortho_file = Path(proj_file).parents[1] / "multispec" / "level1_proc/" / (
                    Path(proj_file).stem + "_" + "multispec_ortho_" + str(res_xy).split('.')[1] + ".tif")

            chunk.exportRaster(path=str(ortho_file), resolution_x=res_xy, resolution_y=res_xy,
                               image_format=Metashape.ImageFormatTIFF,
                               raster_transform=Metashape.RasterTransformValue,
                               source_data=Metashape.OrthomosaicData, image_compression=compression)
            print("Exported orthomosaic: " + str(ortho_file))

        print("Multispec chunk processing complete!")


def main():
    label = "Metashape processing menu"
    Metashape.app.removeMenuItem(label)
    Metashape.app.addMenuItem(label, startProc)
    print("To start the image processing workflow, press {}".format(label))

main()