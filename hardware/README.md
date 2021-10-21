# Hardware module

This module contains all the files that are needed to operate the dedicated read-out device. Additionaly, some files are included that can be used to test if the hardware is working properly. For a full building guide of the device, visit the contribution page of our wiki ([click here](/[modeling/env_requirements](https://2021.igem.org/Team:TUDelft/Contribution))).

## Usage

To run the source code you should install Python and install the required libraries. The required libraries can be found in the folder `hardware/env_requirements` ([click here](/modeling/env_requirements)) and can be installed with Python pip: use `requirements.txt` and install the libraries with the help of *pip*. For help, see [here](https://pip.pypa.io/en/stable/user_guide/#requirements-files).

## File description

A description for files that were used for analysis and text on our wiki page. The filename is shown and it is briefly explained what each file does.

- `file_helpers_test.py` - file to try out the functions that are used to save and read data from csv files.
- `file_helpers.py` - functions to save and read the data of the absorbance and temperature measurements.
- `helpers.py` - helper functions to enable measurements by the read-out device.
- `measurement_main.py` - file to run a temperature-controlled absorbance measurement over time, including an option to pre-heat the test box.
- `measurement_absorbance.py` - file to run an absorbance measurement over time.
- `plots.py` - .
- `temperature_tuning.py` - file to tune the PID-controller.
- `test_i2c_multiplexer.py` - file to try out if the connection and communication with the sensors via the multiplexer functions correctly.
- `test_spi.py` - file to try out the temperature measurement.