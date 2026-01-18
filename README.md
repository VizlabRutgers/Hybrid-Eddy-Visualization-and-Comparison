# Hybrid-Eddy-Visualization-and-Comparison
This is the Eddy visualization code for the hybrid eddy detection and also compare with winding angle and hybrid approach in this paper: https://doi.org/10.48550/arXiv.2305.08229

![plot](./method_comparison.jpg)

![plot](./more_datasets.jpg)

For the comarison, you need the hybrid eddy detection, OW eddy detection and winding angle eddy detection data.

The hybrid eddy detection approach is here: https://github.com/VizlabRutgers/Hybrid-Eddy-detection  
The winding angle detection approach is here: https://github.com/dougcahl/eddy_identification_winding  
The OW eddy detection approach is here: https://github.com/VizlabRutgers/Feature_Tracking  

The dataset presented here in the demo is released here: https://kaust-vislab.github.io/SciVis2020/



# Usage
1. Download the matlab code
2. Create the 3D seabed information (we use Red Sea dataset as a example)
3. Run Multiple eddy detection approaches.
4. Set the data path and configuration at the beginning of the eddyStat.m file.
5. Run the matlab program (It can create video)

See the video for introduction:  
https://drive.google.com/file/d/1t_ST1O7JM3gmUgCUgspKjc9qZZg4cWWJ/view?usp=drive_link
