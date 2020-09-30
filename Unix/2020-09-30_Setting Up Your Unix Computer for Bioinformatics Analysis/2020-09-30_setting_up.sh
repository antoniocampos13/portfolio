# OBJECTIVE: Set up Unix system to use Bioinformatics programs and tools
# SYSTEM: Ubuntu 20.04 LTS (Focal Fossa) 
# Original Post: https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html#setting-up-your-unix-computer-for-bioinformatics-analysis

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

### same effect
./Miniconda3-latest-Linux-x86_64.sh

## Initialize conda:

miniconda3/condabin/conda init

## Close the terminal and open it again. Check it installed properly by invoking miniconda:
conda

## Configuring miniconda channels

conda config --add channels bioconda
conda config --add channels conda-forge

# Create an environment for Bioinformatics programs

conda create -y --name bioenv python=3.8

# Activating an environment
conda activate bioenv

## Installing programs
# Download the `bioenv.txt` file in my GitHub repository. This file contains a selection of most used Bioinformatics programs (hat tip to [Dr. István Albert](https://www.biostarhandbook.com/index.html))

curl 
```bash
cat bioenv.txt | xargs conda install -y
```

## Backing up and restoring your environment configuration

Miniconda has a special command to backup your environment configuration. **Activate** (if needed) the environment you want to backup and enter the command:

```bash
conda env export > bioenv.yml
```

It will result in a `YAML` file in the current working folder containing all configurations in your environment. Again, I named the file `bioenv.yml` but you can choose whatever you like. Note that if you already have a `bioenv.yml` in your directory, it will be overwritten, so be careful.

To restore this environment in your computer, or on other computer, first install miniconda again, and then use the command:

```bash
conda env create -f bioenv.yml
```

The `-f` flag means you are creating an environment using the configurations in the `bioenv.yml` file. The first line of the yml file sets the new environment's name, so you can change it in the file if you like. It will also restore the channels configured in the previous installation of miniconda.

## Conclusion

This is how I configured my system so I could use the major Bioinformatics tools out there. In summary, I:

* Prepared an Unix (Ubuntu) system;
* Installed miniconda, a environment manager;
* Configured channels so I could retrieve desired software;
* Created an environment, showed how to activate and deactivate it, and finally installed software in it;
* Showed how to backup your environment for safekeeping or share with others.

In future posts I will demo some uses of the programs I installed in the new environment.

## References

[Trying the New WSL 2. It's Fast! (Windows Subsystem for Linux) | DigitalOcean](https://www.digitalocean.com/community/posts/trying-the-new-wsl-2-its-fast-windows-subsystem-for-linux)

[Miniconda](https://conda.io/miniconda.html)

[Miniconda installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

[Miniconda & Conda documentation](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

[Managing channels; conda 4.8.4.post65+1a0ab046 documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html)

[The Biostar Handbook: 2nd Edition](https://www.biostarhandbook.com/index.html)