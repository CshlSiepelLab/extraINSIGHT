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

# Covariate directory
COVARIATE_DIR = config["covariates"]["DIRECTORY"]

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
COVARIATE_FILE = os.path.join(COVARIATE_DIR,"unified_covariate_file.bed.gz")

# Set-up wildcards
id = glob_wildcards(os.path.join(HPC_WORKDIR,'jobs','window_chunk.bed.{id}')).id

# Cluster conda environment
cluster_env = os.path.join(HPC_WORKDIR,"environment.yml")

localrules:
    all

rule all:
    input:
        final_mutation_rates = expand(os.path.join(HPC_WORKDIR,'jobs_output','mutation_rates.{id}.bed.gz'), id = id)        
    
# Compute a local adjustment factor for each mb block by calibrating the predicted number of
# mutations from the global model against the observed number of mutations
rule recalibrate_local_mutation_model:
    input:
        global_glm = os.path.join(HPC_WORKDIR,'glm_slim.RData'),
        chrom_size_file = os.path.join(HPC_WORKDIR, CHROM_SIZE_FILE),
        neutral_file = os.path.join(HPC_WORKDIR,NEUTRAL_FILE),
        all_mutations = os.path.join(HPC_WORKDIR,'gnomad_all_potential_mutation.bed.gz'),
        mutation_freq_table = os.path.join(HPC_WORKDIR,'mutation_frequency_table.txt'),
        covariate_file = os.path.join(HPC_WORKDIR, COVARIATE_FILE),
        local_model_script = os.path.join(HPC_WORKDIR,'local_mutation_regression.R'),
        local_windows = os.path.join(HPC_WORKDIR,'jobs','window_chunk.bed.{id}')
    output:
        final_mutation_rates = os.path.join(HPC_WORKDIR,'jobs_output','mutation_rates.{id}.bed.gz'),
        block_status = os.path.join(HPC_WORKDIR,'jobs_output','block_status.{id}.bed')
    shell:
        """
        cd {HPC_WORKDIR}
        # run local regression
        # MIN_N_MUTATION: windows with mutations less that this number are ignored.
        Rscript {input.local_model_script} -b {input.local_windows} -m {input.all_mutations} -n {input.neutral_file} -g {input.global_glm}\
        --flanking-size {FLANKING_SIZE} --min-coverage {MINIMAL_COVERAGE} --max-frequency  {AF_CUTOFF} --min-n-mutation {MIN_N_MUTATION}\
        -k {input.mutation_freq_table} -c {input.covariate_file} -o {output.final_mutation_rates} -l {output.block_status}
        """
