## drone_metashape: RGB and Multispectral imagery processing in Agisoft Metashape Pro
These workflows were developed for RGB and multispectral imagery collected simultaneously on DJI Matrice 300 platform. These scripts are designed for use with RGB imagery acquired with DJI Zenmuse P1 camera and multispectral imagery from MicaSense RedEdge-MX or Dual sensors. For more information please refer to the [Drone Data Collection Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 

One of the following workflows can be used to generate co-registered RGB and multispectral orthomosaics. For a summary of the procesisng steps and information on Agisoft Metashape setup, please refer to the [Drone RGB and Multispectral Processing Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 
1. Automated processing workflow  
metashape_proc.py run from Metashape GUI to generate orthomosaics. User input is required to select images for MicaSense reflectance calibration. 

2. Step-by-step procesing using the Metashape GUI  
metashape_only_upd_cam_pos.py run as a part of the step-by-step workflow to update camera positions in the Metashape project. 

3. Custom widget within Metashape  
This script is only a prototype of an automated processing workflow through a custom widget within Metashape. For this option testing/metashape_proc_widget.py and upd_micasense_pos.py must be copied to C:\Users\<username>\AppData\Local\Agisoft\Metashape Pro\scripts\. For more information refer to Appendix 5 in the [Drone RGB and Multispectral Processing Protocol](https://www.tern.org.au/field-survey-apps-and-protocols/). 


**Funding**: This project was funded by TERN Landscapes  
**Authors**: Poornima Sivanandam, Darren Turner, Arko Lucieer, School of Geography, Planning and Spatial Sciences, University of Tasmania  
**Acknowledgements**: TERN Landscapes, TERN Surveillance, TERN Data Services
