#! /usr/bin/env bash

## Setting up snakemake in this way is crucial so that conda is built and activated correctly, otherwise
## the job dispatch to the cluster will fail as the jobscripts are constructed with a version of python
## that is not accessible from the cluster. --conda-prefix must be a location that is accessible to both
## the cluster and the server!!!

## This makes it so that conda works in a shell session
eval "$(conda shell.bash hook)"
## Activates the extraINSIGHT environment
conda activate extraINSIGHT
## Some configurations
config="../config_files/grch38_neutral_config.yml"
hpc_workdir=`grep "HPC_WORKDIR" ../config_files/grch38_neutral_config.yml | awk '{print $2}'`
conda_dir=${hpc_workdir}/.conda
## Make the conda directory
mkdir -p ${conda_dir}
## Start the pipeline
echo "nohup snakemake --use-conda --conda-prefix ${hpc_workdir}/.conda --configfile=${config} &> logs/grch38.log &"
nohup snakemake --use-conda --conda-prefix ${hpc_workdir}/.conda --configfile=${config} &> logs/grch38.log &
