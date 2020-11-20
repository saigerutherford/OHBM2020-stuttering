# bigdata-stuttering

### Code to support

**Leveraging big data for classification of children who stutter from fluent peers.** 

Saige Rutherford, Mike Angstadt, Chandra Sripada, Soo-Eun Chang. bioRxiv 2020.10.28.359711; doi: https://doi.org/10.1101/2020.10.28.359711

The raw data files necessary for this analysis are restricted (ABCD data). An NDA study (DOI:	10.15154/1520500) has been created for the ABCD data used in this work and can be used to understand which (release 1.1) subjects were included in our analyses. Projected PCA components and other necessary (fully anonymized) variables have been saved in a matlab workspace variable (`data/bigdata-cws.mat`) and can be used only for the purposes of reproducing the results in our [paper](https://www.biorxiv.org/content/10.1101/2020.10.28.359711v1.abstract). 

### Instructions for setup

1. Clone this repository to your computer. `git clone https://github.com/saigerutherford/bigdata-stuttering.git ./bigdata-stuttering`
2. Download files from [Google Drive](https://drive.google.com/drive/folders/1-P3W28lMAmbFbOyHvKYeJ3SX2dOrjOdS?usp=sharing) because they are over GitHub's 100Mb file size limit. Untar the files and add them to the cloned git repo (Hint use: `tar xvf file.tar.gz`).
3. Run `code/s01_PredictiveModel.m` (Note: All analysis was run using MATLAB R2015b and has not been tested in other verisons)
4. Permutations can be reproduced by running `code/s04_Permutations.m`
5. Step away from your computer & do a happy dance or go get a snack. Be creative. 

## Note:

`code/s00_ABCD.m`, `code/s00_HCP.m`, `code/s00_InSample.m` & `code/s03_Keep2Networks.m` are meant for code reviewing purposes only. These scripts contain hard-coded paths to data that are unable to be shared due to restrictive data sharing agreements. Thus these scripts will not be able to be run in a new environment, and are shared to provide transparency for how the data was preprocessed and organized prior to running predictive models. 

This repository is distributed under a [CC-BY-4.0 license](https://creativecommons.org/licenses/by-sa/4.0/).
