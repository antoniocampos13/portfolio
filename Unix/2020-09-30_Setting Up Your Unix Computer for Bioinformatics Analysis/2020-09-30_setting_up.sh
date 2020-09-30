#!/bin/bash

# OBJECTIVE: Set up Unix system to use Bioinformatics programs and tools
# OPERATIONAL SYSTEM: Ubuntu 20.04 LTS (Focal Fossa) 
# ORIGINAL POST: https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html#setting-up-your-unix-computer-for-bioinformatics-analysis

# Updating the system
sudo apt-get update
sudo apt-get upgrade

# Installing useful libraries to be sure that all future libraries will be installed and work properly. Some of these (e.g. default-jdk, the Java libraries), may already be installed in your system, but just to ensure:
sudo apt-get install -y curl unzip build-essential ncurses-dev
sudo apt-get install -y byacc zlib1g-dev python-dev git cmake
sudo apt-get install -y default-jdk ant

# Installing (mini)conda
## Download miniconda installer
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

## Go to the folder where the installer was downloaded and run the script:
bash Miniconda3-latest-Linux-x86_64.sh

## Initialize conda:
miniconda3/condabin/conda init

## Close the terminal and open it again. Check it installed properly by invoking miniconda:
conda

# Configuring miniconda channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Create an environment for Bioinformatics programs
conda create -y --name bioenv python=3.8

# Activating an environment
conda activate bioenv

# Installing Bioinformatics programs
## Download the `bioenv.txt` file in my GitHub repository. This file contains a selection of most used Bioinformatics programs. Hat tip to Dr. István Albert https://www.biostarhandbook.com/index.html
## Otherwise, you also may clone the repository via git.
curl https://raw.githubusercontent.com/antoniocampos13/portfolio/master/Unix/2020-09-30_Setting%20Up%20Your%20Unix%20Computer%20for%20Bioinformatics%20Analysis/bioenv.txt > bioenv.txt

## Find the folder where the file was downloaded and:
cat bioenv.txt | xargs conda install -y

# Backing up and restoring your environment configuration
## Note that if you already have a `bioenv.yml` in your directory, it will be overwritten, so be careful.
## This file is also on the repository already.
conda env export > bioenv.yml

## To restore this environment in your computer, or on other computer, first install miniconda again, and then:
conda env create -f bioenv.yml

# Deactivating an environment
## Remember: you'll have to activate the environment every time you need to use it. 
conda deactivate