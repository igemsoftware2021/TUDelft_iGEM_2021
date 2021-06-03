# Sequence analysis pipeline
Developer: TU Delft iGEM team 2021

Questions to: igem@tudelft.nl

## Setup
### Windows
Download VirtualBox and install a Linux virtual machine. The pipeline was created tested on __Ubuntu 20.04__.
Step by step instructions for the installation of the VirtualBox and Ubuntu: https://itsfoss.com/install-linux-in-virtualbox/.

Follow the steps for Linux on the the VirtualBox.

### MacOS
Same steps as for Linux.

### Linux
The pipeline was created and tested on __Ubuntu 20.04__.

The following steps need to be done on either Linux or MacOS.



The README file for the sequence analysis pipeline.


How to run it
first run "snakemake --cores 1 --use-conda --conda-create-envs-only" This will create the environments needed for the running of the program

Once the environments have been downloaded and installed you can run it using the normal conda commands, but just add the "--use-conda" flag at the end
e.g. : snakemake --cores 2 --use-conda