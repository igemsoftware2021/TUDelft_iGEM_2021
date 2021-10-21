# Modeling

This folder contains all the files used for our model implementation, sensitivity analyses and figure creation. The data is not shown, since the total file size is too large.
## Usage
To run the source code you should install Python and install the required libraries. The required libraries can be found in the folder `modeling/env_requirements` ([click here](/modeling/env_requirements)) and can be installed in two different ways:

1. Python pip: use `requirements.txt` and install the libraries with the help of **pip**. For help, see [here](https://pip.pypa.io/en/stable/user_guide/#requirements-files).
2. Conda package manager: use `environment.yaml` and the conda package manager to create a new environment with the required packages. For help, see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

## File description
A description for files that were used for analysis and text on our wiki page. The filename is shown and what the code in the file does is briefly explained.

- `models.py` - implementation of our model.
- `models_area.py` - implementation of our model for calculating of area between two curves corresponding to two different vitamin concentrations.
- `morris_method.py` - functions to do Morris method on the product concentration as output and the area as output. Also contains functions to write the data and descriptions to files.
- `morris_model_prokaryotic.py` - functions to run the sensitivy analyses with the help of the Morris method.
- `plots.py` - functions to make the figures for the modeling page of the wiki.
- `animations.py` - functions to make the animations for the wiki.
- `plot_helpers.py` - functions that are needed for the plotting.
- `model_timing.py` - file to benchmark model implementations of pure Python + NumPy and Numba.

