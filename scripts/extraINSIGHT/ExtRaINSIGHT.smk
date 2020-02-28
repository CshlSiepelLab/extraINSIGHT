import os

# Configuration parameters
out_dir = config["out_dir"] # directory where output goes
bed = config["bed"] # bed file of query regions
max_sites = config["max_sites"]
mutation_rate = config["mutation_rate"]
exclude_chr = config["exclude_chr"] # this should be a regex expression of chromosomes to filter eg. chrX|chr1

# Other useful parameters
core_filename = os.path.splitext(os.path.basename(bed))[0]

# Reused file names
complete_flag = os.path.join(out_dir, "complete.flag")
subsampled_bed = os.path.join(out_dir, core_filename + "_subsampled.bed")
filtered_bed = os.path.join(out_dir, core_filename + "_filtered.bed")


rule all:
    input:
        complete_flag


# Subsample sites down so that the expected number of sites is equal to the max number of sites
rule subsample_sites:
    input:
        bed = bed
    output:
        sub_bed = subsampled_bed,
        filtered_bed = filtered_bed
    params:
        max_sites = max_sites,
        exclude_chr_regex = exclude_chr
    shell:
        """
        grep -vE {params.exclude_chr_regex} {input.bed} | sort-bed - > {output.filtered_bed}
        ratio=`awk -v max_sites={params.max_sites}  '{{sum = sum + $3 - $2}}END{{print sum}}' {output.filtered_bed}`
        if [ ratio -ge 1 ]; then
            bedops --chop 1 {output.filtered_bed} | awk -v ratio=ratio '{{if (rand() < ratio) print $0}}' | mergeBed -i - > {output.sub_bed}
        else
            sort-bed {output.filtered_bed} | mergeBed -i - > {output.sub_bed}
        fi
        """
# Annotate sites    
rule annotate_sites:
    input:
        sub_bed = subsampled_bed,
        mutation_rate = mutation_rate
    output:
        complete = temp(touch(complete_flag))
    shell:
        """
        tabix {input.mutation_rate} -R {input.sub_bed} | intersectBed -a - -b {input.sub_bed} -sorted | cut -f5,7
        """
