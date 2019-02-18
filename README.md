# Code to reproduce analyses from Luyckx et al. (2019) and experiment files.

Structure:

* Analysis
  - Fig*: reproduces the mentioned figure
  - Preprocessing_pipeline_*: pipelines to preprocess raw EEG data
  - Other auxiliary files
  - Requires:
    - RSA toolbox (http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/)
    - eeglab (https://sccn.ucsd.edu/eeglab/index.php)

* Data
  - Relevant data folders can be downloaded from https://doi.org/10.5061/dryad.7k7s800
  - Contains other empty folders where newly generated data is stored

* Experiments
  - Donkey_EEG: bandit task
  - Numbers_EEG: numerical task
  - Requires:
    - Psychtoolbox-3

* Figures
  - Empty folder that will contain subfolder for generated figures

* Functions
  - cbrewer: extra color maps
  - EEG_preproc_curry: functions to perform preprocessing
  - myfunctions: personally written auxiliary functions
  - rsa: extra functions to perform RSA
  - thirdparty: other functions necessary to run analyses
