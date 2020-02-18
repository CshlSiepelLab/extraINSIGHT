import os.path
import glob
import random
import string


############################################################################
## NEUTRAL MUTATION MODEL ESTIMATION PIPELINE - Cluster subset
## This pipeline handles the subset of the local mutation model correction
## that is run on a cluster so that it can be copied over to the cluster working
## directory and initialize from there. Handles the problem where some networked
## storage may no be accessible from the originating directory of the pipeline
## that triggers this script.
############################################################################

configfile: "configfile.yml"

# Various cutoffs to be used in the pipeline
HIGH_COVERAGE = config["cutoffs"]["HIGH_COVERAGE"] # coverage cutoff for high quality region
MINIMAL_COVERAGE = config["cutoffs"]["MINIMAL_COVERAGE"] # minimal coverage for variant calling
AF_CUTOFF = config["cutoffs"]["AF_CUTOFF"] # maximal allele frequency cutoff for rare alleles
FLANKING_SIZE = config["cutoffs"]["FLANKING_SIZE"] # size of flanking region for local mutation model
MIN_N_MUTATION = config["cutoffs"]["MIN_N_MUTATION"] # minimum number of mutations in a bin

# Output directories
HPC_WORKDIR = config["cluster"]["HPC_WORKDIR"]

# Reference genome associated variables
CHROM_SIZE_FILE = config["reference"]["CHROM_SIZE_FILE"] # chromosome size file
NEUTRAL_FILE = config["reference"]["NEUTRAL_FILE"] # bed file of neutral regions
GC_COVARIATE_FILE = config["covariates"]["GC_COVARIATE_FILE"] # GC covariate file

# Set-up wildcards
id = glob_wildcards(os.path.join(HPC_WORKDIR,'jobs','window_chunk.bed.{id}')).id

# Cluster conda environment
cluster_env = os.path.join(HPC_WORKDIR,"environment.yml")

localrules:
    all

rule all:
    input:
        final_mutation_rates = expand(os.path.join(HPC_WORKDIR,'jobs_output','mutation_rates.{id}.bed.gz'), id = id)        
    output:
        touch(os.path.join(HPC_WORKDIR,'proc_complete.txt'))
 
## Pre-initialize conda environment so that the parallel jobs don't try to do it many times
rule initialize_conda_env:
    conda:
        cluster_env
    output:
        flag = touch(os.path.join(HPC_WORKDIR,"conda_env_created.flag"))
    shell:
        """
        cd {HPC_WORKDIR}
        echo "setting up conda environment for cluster jobs..."
        """
        
## Need to preinstall the conda environment on the cluster for this to work, otherwise just make sure
## the dependencies are installed for now, will try and get it to work later        
rule recalibrate_local_mutation_model:
    input:
        flag = rules.initialize_conda_env.output.flag,
        global_glm = os.path.join(HPC_WORKDIR,'glm_slim.RData'),
        chrom_size_file = os.path.join(HPC_WORKDIR, CHROM_SIZE_FILE),
        neutral_file = os.path.join(HPC_WORKDIR,NEUTRAL_FILE),
        all_mutations = os.path.join(HPC_WORKDIR,'gnomad_all_potential_mutation.bed.gz'),
        mutation_freq_table = os.path.join(HPC_WORKDIR,'mutation_frequency_table.txt'),
        covariate_file = os.path.join(HPC_WORKDIR,GC_COVARIATE_FILE),
        local_model_script = os.path.join(HPC_WORKDIR,'local_mutation_regression.R'),
        par_model_script = os.path.join(HPC_WORKDIR,'build_local_mutation_model_parallel.pl'),
        local_windows = os.path.join(HPC_WORKDIR,'jobs','window_chunk.bed.{id}')
    output:
        final_mutation_rates = os.path.join(HPC_WORKDIR,'jobs_output','mutation_rates.{id}.bed.gz')
    conda:
        cluster_env
    shell:
        """
        cd {HPC_WORKDIR}
        # run local regression
        # MIN_N_MUTATION: windows with mutations less that this number are ignored.
        perl {input.par_model_script} {input.local_windows} {input.all_mutations} {input.neutral_file} {input.global_glm}\
        {FLANKING_SIZE} {MINIMAL_COVERAGE} {AF_CUTOFF} {MIN_N_MUTATION}\
        {input.mutation_freq_table} {input.covariate_file} | bgzip -c > {output.final_mutation_rates}
        """
