# Sequence analysis pipeline
Developer: TU Delft iGEM team 2021

Questions to: igem@tudelft.nl

## Setup
### Windows
Download VirtualBox and install a Linux virtual machine. The pipeline was created tested on __Ubuntu 20.04__.
Step by step instructions for the installation of the VirtualBox and Ubuntu: https://itsfoss.com/install-linux-in-virtualbox/.

Now follow the steps for __Linux/MacOS__ in the VirtualBox.

### Linux/MacOS
The pipeline was created and tested on __Ubuntu 20.04__. Here are the steps one need to take:
1. Install Miniconda or Anaconda with Python version 3.9.* or higher. If you have not much experience with Anaconda or Python we recommend downloading Anaconda.
2. Open the terminal.
    * Linux: `CTRL + ALT + t`
    * MacOS: Press `Cmd + Space` to open Spotlight search, and type terminal and hit return.
3. Write `conda activate` to activate the conda environment. (Link to getting started with conda: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html.)
4. 

The following steps need to be done on either Linux or MacOS.



The README file for the sequence analysis pipeline.


How to run it
first run "snakemake --cores 1 --use-conda --conda-create-envs-only" This will create the environments needed for the running of the program

Once the environments have been downloaded and installed you can run it using the normal conda commands, but just add the "--use-conda" flag at the end
e.g. : snakemake --cores 2 --use-conda