# tissue_respiration
This repository contains data and analysis scripts corresponding to the article
**Short-term thermal acclimation of dark respiration is greater in non-photosynthetic 
than in photosynthetic tissues** by Smith et al. in *Annals of Botany - Plants*.

## Description of repository folders

### Raw Data
The [raw_data](raw_data) folder contains raw gas exchange data for respiration of 
[leaves](raw_data/leaf_raw.csv), [stems](raw_data/stem_raw.csv), and [roots](raw_data/root_raw.csv)
as well as [photosynthetic and leaf trait data](gc_data_merged.csv). Please refer to
[https://github.com/SmithEcophysLab/PU_GrowthChamber](https://github.com/SmithEcophysLab/PU_GrowthChamber)
for more information about the photosynthetic data.

### Curve Fitting Script
The [curve_fitting_script](curve_fitting_script) folder contains an 
[R script](curve_fitting_script/stem_root_curvefitting.R) to fit the instantaneous respiration
temperature response curves from the [raw data](raw_data).

### Instantaneous Temperature Response Curve Fits
The [tresp_curve_fits](tresp_curve_fits) folder contains the fits for each instantaneous
temperature response of [leaf](tresp_curve_fits/leaf_fits_v3.csv), [stem](tresp_curve_fits/stem_fits_v3.csv),
and [root](tresp_curve_fits/root_fits_v3.csv) respiration from the curve fitting
[R script](curve_fitting_script/stem_root_curvefitting.R).

### Analysis Script
The [analysis_script](analysis_script) folder contains an 
[R script](analysis_script/tissue_respiration_analysis.R) to run the analyses presented in the
Smith et al. article.