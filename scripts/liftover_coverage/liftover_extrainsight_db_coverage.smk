import os.path
import glob
import random
import string

## Retrieve variables from configfile
ROOT_DIR = config["rootdir"]
OUTPUT_DIR = os.path.join(ROOT_DIR, config["output"]["MUTMOD_OUTPUT_DIR"])
GENOME_NAME = config["reference"]["GENOME_NAME"] # name of the genome

## Set some important variables
CHAIN_DIR = "../../data/grch38/liftover/"
if GENOME_NAME == "grch38":
    CHAIN_FILE = os.path.join(CHAIN_DIR, "hg38ToHg19.over.chain.gz")
    TARGET_GENOME = "grch37"
if GENOME_NAME == "grch37":
    CHAIN_FILE = os.path.join(CHAIN_DIR, "hg19ToHg38.over.chain.gz")
    TARGET_GENOME = "grch38"

## Conda environment
local_env = os.path.join("..","..","environment.yml")

def aggregate_success(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the create_local_windows step
    '''
    # get the 'job_directory' variable from the checkpoint
    checkpoint_output_dir = checkpoints.split_coverage.get(**wildcards).output.out_dir
    # Get suffixes of split output files
    id = glob_wildcards(os.path.join(OUTPUT_DIR,'split_coverage','{id}.bed')).id
    # sort by wildcards as files are internally sorted and the split is ordered so ordering the wildcards
    # preserves the order for concatination
    # Generate the wildcards from the output of the checkpoint rule
    return expand(os.path.join(OUTPUT_DIR,'split_liftover',
                               f'{{id}}.{GENOME_NAME}.to.{TARGET_GENOME}.success.bed'),
                  id = id)

def aggregate_failure(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the create_local_windows step
    '''
    # get the 'job_directory' variable from the checkpoint
    checkpoint_output_dir = checkpoints.split_coverage.get(**wildcards).output.out_dir
    # Get suffixes of split output files
    id = glob_wildcards(os.path.join(OUTPUT_DIR,'split_coverage','{id}.bed')).id
    # sort by wildcards as files are internally sorted and the split is ordered so ordering the wildcards
    # preserves the order for concatination
    # Generate the wildcards from the output of the checkpoint rule
    return expand(os.path.join(OUTPUT_DIR,'split_liftover',
                               f'{{id}}.{GENOME_NAME}.to.{TARGET_GENOME}.fail.bed'),
                  id = id)

rule all:
    input:
        os.path.join(OUTPUT_DIR, f"{TARGET_GENOME}.syntenic.site.coverage.bed.gz"),
        os.path.join(OUTPUT_DIR, f"{GENOME_NAME}.syntenic.site.coverage.bed.gz")

checkpoint split_coverage:
    input:
        source_coverage = os.path.join(OUTPUT_DIR,"final_mutation_site_coverage.bed.gz")
    output:
        out_dir = temp(directory(os.path.join(OUTPUT_DIR, "split_coverage"))),
        lo_dir = temp(directory(os.path.join(OUTPUT_DIR, "split_liftover")))
    conda:
        local_env
    shell:
        """
        mkdir -p {output.out_dir}
        mkdir -p {output.lo_dir}
        zcat {input.source_coverage} |\
        awk 'BEGIN{{OFS = "\t"}}{{print > "{output.out_dir}/"$1".bed"}}'
        """
        
rule liftover_coverage:
    input:
        os.path.join(OUTPUT_DIR, "split_coverage"),
        os.path.join(OUTPUT_DIR, "split_liftover"),
        db = os.path.join(OUTPUT_DIR, "split_coverage","{id}.bed")
    output:
        success = os.path.join(OUTPUT_DIR, "split_liftover",
                               f"{{id}}.{GENOME_NAME}.to.{TARGET_GENOME}.success.bed"),
        fail = os.path.join(OUTPUT_DIR, "split_liftover",
                            f"{{id}}.{GENOME_NAME}.to.{TARGET_GENOME}.fail.bed")
    conda:
        local_env
    shell:
        """
        bedops -w 1 {input.db} | liftOver /dev/stdin {CHAIN_FILE} {output.success} {output.fail}
        bedops -m {output.success} > {output.success}.tmp
        mv {output.success}.tmp {output.success}
        bedops -m {output.fail} > {output.fail}.tmp
        mv {output.fail}.tmp {output.fail}
        """

rule concatonate_valid_target_genome_sites:
    input:
        os.path.join(OUTPUT_DIR, "split_liftover"),
        gathered_successes = aggregate_success,
    output:
        valid_target = os.path.join(OUTPUT_DIR, f"{TARGET_GENOME}.syntenic.site.coverage.bed.gz")
    conda:
        local_env
    shell:
        """
        cat {input.gathered_successes} | sort-bed - | bedops -m - | gzip -c - > {output.valid_target} 
        """
        
## Get valid sites on source genome by removing the sites that failed to lift over
rule concatonate_valid_source_genome_sites:
    input:
        os.path.join(OUTPUT_DIR, "split_liftover"),
        source_coverage = os.path.join(OUTPUT_DIR,"final_mutation_site_coverage.bed.gz"),
        gathered_failure = aggregate_failure
    output:
        valid_source = os.path.join(OUTPUT_DIR, f"{GENOME_NAME}.syntenic.site.coverage.bed.gz")
    conda:
        local_env
    shell:
        """
        set +o pipefail;
        cat {input.gathered_failure} | grep -v "#" | sort-bed - | bedops -m - |\
        bedops -d <(zcat {input.source_coverage}) - | gzip -c -  > {output.valid_source} 
        """

