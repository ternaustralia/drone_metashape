## drone_metashape: RGB and Multispectral imagery processing in Agisoft Metashape Pro
These workflows were developed for RGB and multispectral imagery collected simultaneously on DJI Matrice 300 platform. These scripts are designed for use with RGB imagery acquired with DJI Zenmuse P1 camera and multispectral imagery from MicaSense RedEdge-MX or Dual sensors. For more information please refer to the [Drone Data Collection Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 

One of the following workflows can be used to generate co-registered RGB and multispectral orthomosaics. For a summary of the processing steps and information on Agisoft Metashape setup, please refer to the [Drone RGB and Multispectral Processing Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 
1. Automated processing workflow  
Uses scripts metashape_proc.py and upd_micasense_pos.py. 
metashape_proc.py run from Metashape GUI to generate orthomosaics. User input is required to select images for MicaSense reflectance calibration. 

2. Step-by-step procesing using the Metashape GUI  
Uses scripts metashape_only_upd_cam_pos.py and upd_micasense_pos.py. 
metashape_only_upd_cam_pos.py run as a part of the step-by-step workflow to update camera positions in the Metashape project. 

3. Prototype only: custom widget within Metashape  
This script is a prototype of an automated processing workflow through a custom widget within Metashape. For this option testing/metashape_proc_widget.py and upd_micasense_pos.py must be copied to C:\Users\<username>\AppData\Local\Agisoft\Metashape Pro\scripts\. For more information please refer to Appendix 5 in the [Drone RGB and Multispectral Processing Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 

4. Other examples<br>
examples\metashape_blockshift.py - code to only perform the blockshift of images in a Metashape chunk using AUSPOS results.  <br>
examples\metashape_proc_p1.py - only process RGB images captured using Zenmuse P1 on **gimbal 1 of dual mount**. **Remove the GPS/INS offset code if P1 was on single mount gimbal.** <br>

**Funding**: This project was funded by TERN Landscapes  
**Authors**: Poornima Sivanandam, Darren Turner, Arko Lucieer, School of Geography, Planning and Spatial Sciences, University of Tasmania  
**Acknowledgements**: TERN Landscapes, TERN Surveillance, TERN Data Services
