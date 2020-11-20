# bigdata-stuttering

### Code to support

Leveraging big data for classification of children who stutter from fluent peers. Saige Rutherford, Mike Angstadt, Chandra Sripada, Soo-Eun Chang. bioRxiv 2020.10.28.359711; doi: https://doi.org/10.1101/2020.10.28.359711

The raw data files necessary for this analysis are restricted (ABCD data). An NDA study (DOI:	10.15154/1520500) has been created for the ABCD data used in this work. However, all data has been anonymized and a shared matlab workspace variable (`data/bigdata-cws.mat`) can be used to re-run the predictive modeling portion. 

### Instructions for setup

1. Clone this repository to your computer. `git clone https://github.com/saigerutherford/bigdata-stuttering.git ./bigdata-stuttering`
2. Run `code/s02_PredictiveModel.m` (Note: All analysis was run using MATLAB R2015b and has not been tested in other verisons)

## Note:

`code/s00_ABCD.m`, `code/s00_HCP.m`, & `code/s00_InSample.m` are meant for code reviewing purposes only. These scripts contain hard-coded paths to data that is unable to be shared due to restrictive data sharing agreements. Thus these scripts will not be able to be run in a new environment, and are shared to provide transparency for how the data was preprocessed and organized prior to running predictive models. 

This repository is distributed under a CC-BY-4.0 license [https://creativecommons.org/licenses/by-sa/4.0/]
