# Hardware module

This module contains all the files that are needed to operate the dedicated read-out device. Additionaly, some files are included that can be used to test if the hardware is working properly. For a full building guide of the device, visit the contribution page of our wiki ([click here](/[modeling/env_requirements](https://2021.igem.org/Team:TUDelft/Contribution))).

## Usage
To run the source code you should install Python and install the required libraries. The required libraries can be found in the folder `hardware/env_requirements` ([click here](/modeling/env_requirements)) and can be installed in two different ways:
1. Python pip: use `requirements.txt` and install the libraries with the help of *pip*. For help, see [here](https://pip.pypa.io/en/stable/user_guide/#requirements-files).
2. Conda package manager: use `environment.yaml` and use the *conda* package manager to create a new environment with the required libraries. For help, see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).


# Modeling

This folder contains all the files used for our model implementation, sensitivity analyses and figure creation. The data is not shown, since the total size of all the files combined is too large to store on GitHub.
## Usage
To run the source code you should install Python and install the required libraries. The required libraries can be found in the folder `modeling/env_requirements` ([click here](/modeling/env_requirements)) and can be installed in two different ways:

1. Python pip: use `requirements.txt` and install the libraries with the help of *pip*. For help, see [here](https://pip.pypa.io/en/stable/user_guide/#requirements-files).
2. Conda package manager: use `environment.yaml` and use the *conda* package manager to create a new environment with the required libraries. For help, see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

## File description
A description for files that were used for analysis and text on our wiki page. The filename is shown and it is briefly explained what each file does.

- `file_helpers_test.py` - file to try out the functions that are used to save and read data from csv files.
- `file_helpers.py` - functions to 
- `helpers.py` - a
- `measurement_absorbance.py` - a 
- `measurement.py` - a 
- `plots.py`
- `temperature_tuning.py`
- `test_i2c_multiplexer.py`
- `test_spi.py`

- `models.py` - implementation of our model.
- `models_area.py` - implementation of our model for calculating of area between two curves corresponding to two different vitamin concentrations.
- `morris_method.py` - functions to do Morris method on the product concentration as output and the area as output. Also contains functions to write the data and descriptions to files.
- `morris_model_prokaryotic.py` - functions to run the sensitivy analyses with the help of the Morris method.
- `plots.py` - functions to make the figures for the modeling page of the wiki.
- `animations.py` - functions to make the animations for the wiki.
- `plot_helpers.py` - functions that are needed for the plotting.
- `model_timing.py` - file to benchmark model implementations of pure Python + NumPy and Numba.

