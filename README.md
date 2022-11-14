

## Getting started with this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installaed Miniconda3 version `py39_4.12.0` and mamba version `0.15.3`.

```
mamba env create -n tx_gx --file environment.yml
conda activate tx_gx
```

To start the pipeline, run:
```
snakemake --use-conda -j 2
```

### Running this repository on AWS

This repository was executed on an AWS EC2 instance (Ubuntu 22.04 LTS ami-085284d24fe829cd0, t2.2xlarge, 200 GiB EBS `gp2` root storage).
The instance was configured using the following commands:

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba # install mamba for faster software installation.
```

## Papers that inspired this project

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1873-2
https://www.biorxiv.org/content/10.1101/2022.09.25.509396v1
