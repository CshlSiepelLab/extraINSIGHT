#! /usr/bin/env python

import argparse
import os
import shlex
import shutil

# Set textwap width for parser
os.environ['COLUMNS'] = "90"

root_dir = os.path.expanduser("~/Projects/extraINSIGHT")
print(root_dir)

# Specify parser
parser = argparse.ArgumentParser()
parser.add_argument("-b","--bed", dest="bed", type=str,
                    required=True, help="bed file containing query regions")
parser.add_argument("-g","--bed_genome", dest="bed_genome", type=str,
                    required=True,
                    choices=['hg19', 'hg38'],
                    help="the genome version the bed file is from. Either hg19 or hg38")
parser.add_argument("--insight-genome", dest="insight_genome", type=str, default="hg19",
                    choices=['hg19'],
                    help="The genome version that INSIGHT is using. Do not touch unless you are sure.")
parser.add_argument("--extrainsight-genome", dest="extra_insight_genome", type=str, default="hg38",
                    choices=['hg19', 'hg38'],
                    help="The genome version that ExtRaINSIGHT is using. Do not touch unless you are sure.")
parser.add_argument("-o","--out-dir", dest="out_dir", type=str,
                    default = ".", help="results output directory (default = '.')")
parser.add_argument("--exclude-chr", dest="exclude_chr", nargs="+", type=str, default = None,
                    required=False, help="chromosomes to exclude from the analysis")
args = parser.parse_args()

def make_directory(d):
    if not os.path.exists(d):
        os.makedirs(d)
        print("Directory " , d ,  " created ")
    else:    
        print("Directory " , d ,  " already exists")

# Function for calling liftOver to liftover bed files one base at a time, then sorting and merging
# the output. It also creates a bed file of the sites from the original file that were successfully
# mapped over. The function returns the file names for the files containing the intervals that
# successfully mapped over from the perspective of the source and target genome respectively
def sitewise_liftover(bed, chain, out_dir, source_genome, target_genome):
    # Core id for file name generation
    core_name = os.path.splitext(os.path.basename(bed))[0]
    # File of lifted over intervals
    out_mapped = os.path.join(out_dir, core_name + "." + target_genome + "mapped.bed")
    # File of intervals that failed to lift over
    out_unmapped = os.path.join(out_dir, core_name + "." + source_genome + "unmapped.bed")
    # File of intervals that lifted over but on the source genome
    out_mapped_source = os.path.join(out_dir, core_name + "." + source_genome + "mapped.bed")
    # This command chops up input bed file into single base intervals then lifts them over
    cmd_lift = f"bedops --chop 1 {bed} | liftOver /dev/stdin {chain} {out_mapped}.tmp {out_unmapped}"
    # This command sorts and merges the lifted over bases
    cmd_merge = f"sort-bed {out_mapped}.tmp | bedops -m - > {out_mapped}"
    # This command gets the set difference between the input file and the elements that failed to
    # liftOver to create a file of the intervals that successfully lifted over
    cmd_success = f"grep -v '#' {out_unmapped} | bedops -d {bed} - > {out_mapped_source}"
    os.system(cmd_lift)
    os.system(cmd_merge)
    os.system(cmd_success)
    return out_mapped_source, out_mapped

# Create output directory
make_directory(args.out_dir)
     
# Locations of the chain files
hg19_to_hg38_chain = os.path.join(root_dir, "data/grch38/liftover/hg19ToHg38.over.chain.gz")
hg38_to_hg19_chain = os.path.join(root_dir, "data/grch38/liftover/hg38ToHg19.over.chain.gz")

# Make dictionary of chains, first level is source, second level is target
chain_dict = { 'hg19': {'hg38': hg19_to_hg38_chain},
               'hg38': {'hg19': hg38_to_hg19_chain}}

# If the insight and extraINSIGHT genome versions do not match then some lifting over will have to be done
bed_dir = os.path.join(args.out_dir, "bed")
# Create directory to store bed files
make_directory(bed_dir)

# Copy bed file to bed directory and ensure sorted
# Remove excluded chromosomes from file if specified
bed_local = os.path.join(bed_dir, os.path.basename(args.bed))
if args.exclude_chr is not None:
    excl = "|".join(args.exclude_chr)
    os.system(f"awk '$1 !~ /{excl}/' {args.bed} | sort-bed - > {bed_local}")
else :
    os.system(f"sort-bed {args.bed} > {bed_local}")
    
# Set defaults for bed files to use in INSIGHT/ExtRaINSIGHT
i_bed = bed_local
ei_bed = bed_local

# Liftover the bed file to the insight genome if needed
if not args.bed_genome == args.insight_genome:
    chosen_chain = chain_dict[args.bed_genome][args.insight_genome]
    ei_bed, i_bed = sitewise_liftover(bed_local, chosen_chain, bed_dir, args.bed_genome, args.insight_genome)

# Liftover the bed file to the ExtRaINSIGHT genome if needed
if not args.bed_genome == args.extra_insight_genome:
    chosen_chain = chain_dict[args.bed_genome][args.extra_insight_genome]
    i_bed, ei_bed = sitewise_liftover(bed_local, chosen_chain, bed_dir, args.bed_genome, args.extra_insight_genome)

i_bed = os.path.abspath(i_bed)
ei_bed = os.path.abspath(ei_bed)

# Setup insight analysis directories
extra_dir = os.path.abspath(os.path.join(args.out_dir, "ExtRaINSIGHT"))
make_directory(extra_dir)
insight_dir = os.path.abspath(os.path.join(args.out_dir, "INSIGHT"))
make_directory(insight_dir)

cwd = os.getcwd()
# Run ExtRaINSIGHT
os.chdir('../extraINSIGHT/')
cmd_extra = f"""./ExtRaINSIGHT.py -b {i_bed} -o {extra_dir} \
-r /local/projects/extraINSIGHT/results/grch38/gnomad_v3.0/mutation_model/final_mutation_rates.bed.gz"""
os.system(cmd_extra)

os.chdir(cwd)

# Run INSIGHT
cmd_insight = f"../INSIGHT2/INSIGHT2.py -b {i_bed} -o {insight_dir}"
os.system(cmd_insight)
