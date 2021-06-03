# Sequence analysis pipeline
Developer: TU Delft iGEM team 2021
Questions to: igem@tudelft.nl

## Setup
If on *Windows*:
Download VirtualBox and install a Linux virtual machine. The pipeline was created tested on *Ubuntu 20.04*.
Step by step instructions: https://itsfoss.com/install-linux-in-virtualbox/.

The following steps are for Linux and MacOS.



The README file for the sequence analysis pipeline.


How to run it
first run "snakemake --cores 1 --use-conda --conda-create-envs-only" This will create the environments needed for the running of the program

Once the environments have been downloaded and installed you can run it using the normal conda commands, but just add the "--use-conda" flag at the end
e.g. : snakemake --cores 2 --use-conda