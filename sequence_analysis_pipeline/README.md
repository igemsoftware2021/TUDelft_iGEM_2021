[![Snakemake](https://img.shields.io/badge/snakemake-≥6.4.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
# Sequence analysis pipeline
Developer: TU Delft iGEM team 2021

Questions to: igem@tudelft.nl

## Setup
### Windows
Download VirtualBox and install a Linux virtual machine. The pipeline was created and tested on __Ubuntu 20.04__.
Step by step instructions for the installation of the VirtualBox and Ubuntu: https://itsfoss.com/install-linux-in-virtualbox/.

Once you have installed the Linux virtual machine follow the steps for __Linux/MacOS__ in the VirtualBox.

### Linux/MacOS
The pipeline was created and tested on __Ubuntu 20.04__. Here are the steps one need to take:
1. Install Miniconda or Anaconda with Python version 3.9.* or higher. If you have not much experience with Anaconda or Python we recommend installing Anaconda.
2. Open the terminal.
    * Linux: `CTRL + ALT + t`
    * MacOS: Press `Cmd + Space` to open Spotlight search, and type terminal and hit return.
3. Write `conda activate` to activate the conda environment. (Link to getting started with conda: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html.)
4. Now follow the installation steps from snakemake, make sure to do the __full__ installation. The link to the installation steps: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html.
5. Congratulations you now have succesfully installed snakemake! We recommand that you follow the short snakemake tutorial on the site of snakemake. This will help you understand how snakemake works and how our pipeline is built. Here is the link to the tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html.

## Usage
### Step 1: Install workflow
If you simply want to use this workflow, download and extract the latest release. If you intend to modify and further develop this workflow, fork this repository.
### Step 2: Configure the workflow
Configure the workflow to your needs by editing the `config.yaml` and the sample 

### Prepare the data
1. Create a folder called `NGS` in the `data` folder in the `sequence_analysis_pipeline` folder.
2. Insert the forward and reverse Illumina reads in this folder. They will look something like `..._R1_001.fastq.gz` and `..._R2_001.fastq.gz`. File names are important in the pipeline, so you need to rename them to make the pipeline work correctly. At the moment the following patterns are used to retrieve information from the filenames:
    * S... -> selection form which the sequences are from
    * D... -> what was the DRIVER round at the moment of sequencing

    Here is an example, the files `S1_D63_R1_001.fastq.gz` and `S1_D63_R2_001.fastq.gz` indicate that the read sequences are from selection __1__, and the sequencing was done after round __63__ . You can choose to use the same naming convention or change the matching patterns in the `config.yaml` file in the `config` folder.
3. Adjust the file called `ngs_references.csv` in the `data` folder. Change the reference sequences to your own reference sequences. The columns __name__ and __sequence__ will be read in during the running of the pipeline.
4. Go into the `config.yaml` file in the `config` folder and set the name of the dataset. Set the prefixes and suffixes, so whether a prefix or suffix indicates cleaving or uncleaving and what the corresponding name and sequence are of this prefix/suffix. An example config-file can be seen in the `config` folder called `config_example.yaml`.

### Run the software
1. Open the `sequence_analysis_pipeline` folder.
2. Activate your conda snakemake environment via
```
conda activate *yoursnakemake_environmentname*
```
Test your configuration by performing a dry run via
```
snakemake -n
```
_Optional_: The first time you will run the pipeline it is smart to create the needed conda environments beforehand. Let snakemake create the needed environments via
```
snakemake --cores 1 --use-conda --conda-create-envs-only
```
Execute the workflow locally via
```
snakemake --cores $N --use-conda
```
using $N cores.

3. Before we can run the pipeline, snakemake needs to create the needed conda environments. To do this run the following command in the terminal: `snakemake --cores 1 --use-conda --conda-create-envs-only`. Now wait for the environments to be downloaded.
4. Now run the pipeline by using the following command `snakemake --cores $N --use-conda`, where __N__ is the number of cores you let the pipeline use.

## References
Literature used for our sequence analysis pipeline.

* Townshend, B., Xiang, J. Manzanarez, G., Hayden, E. & Smolke, C. *A multiplexed, automated evolution pipeline enables scalable discovery and characterization of biosensors*. http://biorxiv/lookup/doi/10.1101/2020.05.29.117960 (2020) doi:10.1101/2020.05.29.117960
* Zhang, J., Kobert, K., Flouri T. & Stamatakis, A. *PEAR: a fast and accurate Illumina Paired-End reAd mergeR.* Bioinformatics (Oxford, England) 30, 614-20 (2014).
* Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. *Sustainable data analysis with Snakemake.* F1000Res 10, 33.
