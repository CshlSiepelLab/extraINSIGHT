#! /usr/bin/env python

import argparse
import os
import shlex
import shutil
import string
import random
import warnings
import gzip
import io
import subprocess, signal
        
# Set textwap width for parser
os.environ['COLUMNS'] = "90"
root_dir = os.path.expanduser("~/Projects/extraINSIGHT")
#print(root_dir)

# Specify parser
parser = argparse.ArgumentParser()
parser.add_argument("-b","--bed", dest="bed", type = str,
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
parser.add_argument("--max-sites", dest="max_sites", type=int, default = 1e5,
                    required=False, help="maximum number of sites to analyze, if exceeded it will subsample")
args = parser.parse_args()

# Resolve args.out_dir to absolute path
args.out_dir = os.path.abspath(args.out_dir)

def bash_call(cmd, verbose = False):
    cmd_final = "/usr/bin/env bash -c " + f"\"{cmd}\"" + ";"
    if verbose:
        print(cmd_final)
    os.system(cmd_final)

def make_directory(d):
    if not os.path.exists(d):
        os.makedirs(d)
        print("Directory " , d ,  " created ")
    else:    
        print("Directory " , d ,  " already exists")

def random_string(stringLength=8):
    letters = string.ascii_lowercase
    return(''.join(random.choice(letters) for i in range(stringLength)))

def lines(file):
    count = 0
    ext = os.path.splitext(os.path.basename(file))[-1]
    if ext == ".gz":
        gz = gzip.open(file, 'rb')
        f = io.BufferedReader(gz)
        for line in f:
            count += 1
        gz.close()
    else:
        with open(file, 'rt') as f:
            for line in f:
                count += 1        
    return(count)

def count_sites(bed):
    sites = 0
    ext = os.path.splitext(os.path.basename(bed))[-1]
    if ext == ".gz":
        gz = gzip.open(bed, 'rb')
        b = io.BufferedReader(gz)
        ln = b.readline()
        txt = ln.decode("utf8").strip().split("\t")        
        while len(txt) >= 3:
            sites += int(txt[2]) - int(txt[1])
            ln = b.readline()
            txt = ln.decode("utf8").strip().split("\t") 
        gz.close()
    else:
        with open(bed) as b:
            ln = b.readline()
            txt = ln.strip().split("\t")
            while len(txt) >= 3:
                sites += int(txt[2]) - int(txt[1])
                ln = b.readline()
                txt = ln.strip().split("\t")
    return(sites)

def bash_file_handler(file):
    ext = os.path.splitext(os.path.basename(file))[-1]
    if ext == ".gz":
        handler = f"<(zcat {file})"
    else:
        handler = f"<(cat {file})"
    return(handler)

# Function for calling liftOver to liftover bed files one base at a time, then sorting and merging
# the output. It also creates a bed file of the sites from the original file that were successfully
# mapped over. The function returns the file names for the files containing the intervals that
# successfully mapped over from the perspective of the source and target genome respectively
def sitewise_liftover(bed, chain, out_dir, source_genome, target_genome):
    if source_genome == target_genome:
            warnings.warn("source and target genome are equal")
    # Core id for file name generation
    core_name = os.path.splitext(os.path.basename(bed))[0]
    ext = os.path.splitext(os.path.basename(bed))[-1]
    # File of lifted over intervals
    out_mapped = os.path.join(out_dir, core_name + "." + target_genome + ".mapped.bed")
    # File of intervals that failed to lift over
    out_unmapped = os.path.join(out_dir, core_name + "." + source_genome + ".unmapped.bed")
    # File of intervals that lifted over but on the source genome
    out_mapped_source = os.path.join(out_dir, core_name + "." + source_genome + ".mapped.bed")
    if ext == ".gz":
        bed = f"<(zcat {bed})"
    # This command chops up input bed file into single base intervals then lifts them over
    cmd_lift = f"sort-bed {bed} | bedops -w 1 - | liftOver /dev/stdin {chain} {out_mapped}.tmp {out_unmapped}"
    # This command sorts and merges the lifted over bases
    cmd_merge = f"sort-bed {out_mapped}.tmp | bedops -m - > {out_mapped}"
    # This command gets the set difference between the input file and the elements that failed to
    # liftOver to create a file of the source intervals that successfully lifted over
    bash_call(cmd_lift)
    bash_call(cmd_merge)
    if lines(out_unmapped) > 0:
        cmd_success = f"grep -v '#' {out_unmapped} | sort-bed - | bedops -d {bed} - | bedops -m - > {out_mapped_source}"
    else:
        cmd_success = f"sort-bed {bed} | bedops -m - > {out_mapped_source}"
    bash_call(cmd_success)
    os.remove(f"{out_mapped}.tmp")
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
ext = os.path.splitext(os.path.basename(args.bed))[-1]
bed_handler = bash_file_handler(args.bed)
if ext == ".gz":
    bed_local =  os.path.join(bed_dir, os.path.basename(args.bed))
else:
    bed_local = os.path.join(bed_dir, os.path.basename(args.bed) + ".gz")

if args.exclude_chr is not None:
    excl = "|".join(args.exclude_chr)
    bash_call("awk '\$1 !~" + f"/{excl}/' {bed_handler} | sort-bed - | bedops -m - | gzip -c > {bed_local}")
else:
    bash_call(f"sort-bed {bed_handler} | bedops -m - | gzip -c  > {bed_local}")
    
# Set the mutation rate file for ExtRaINSIGHT
if args.extra_insight_genome == "hg38":
    ei_mut_rate_path = os.path.join(root_dir, "results/grch38/gnomad_v3.0/mutation_model/final_mutation_rates.bed.gz")
    ei_mut_cover_path = os.path.join(root_dir, "results/grch38/gnomad_v3.0/mutation_model/final_mutation_site_coverage.bed.gz")
elif args.extra_insight_genome == "hg19":
    ei_mut_rate_path = os.path.join(root_dir, "results/grch37/gnomad_v2.1/mutation_model/final_mutation_rates.bed.gz")
    ei_mut_cover_path = os.path.join(root_dir, "results/grch37/gnomad_v2.1/mutation_model/final_mutation_site_coverage.bed.gz")

## To ensure that the same exact sites are being analyzed by INSIGHT2 and ExtRaINSIGHT:
# 1) Ensure bed file in co-ordinate system of ExtRaINSIGHT
# 2) Filter out sites not in ExtRaINSIGHT database
# 3) LiftOver to INSIGHT co-ordinate system
# 4) Keep the sites that lifted over from the ExtRaINSIGHT co-ordinate system for analysis

initial_sites = count_sites(bed_local)
print(f"Sites in intial bed file: {initial_sites}")
# if greater than max sites, subsample
if args.max_sites is not None and initial_sites > args.max_sites:
    print("subsampling...")
    fraction = float(args.max_sites)/initial_sites
    tmp_out = os.path.join(bed_dir, os.path.basename(args.bed)) + ".tmp"
    bed_handler = bash_file_handler(bed_local)
    subsample_cmd = f"bedops --chop 1 {bed_handler} | awk -v frac={fraction} '{{if(rand() <= frac) print \$0}}' | gzip -c > {tmp_out}"
    bash_call(subsample_cmd, True)
    shutil.move(tmp_out, bed_local)

subsampled_sites = count_sites(bed_local)
print(f"Sites after subsampling: {subsampled_sites}")
    
# Liftover the bed file to the ExtRaINSIGHT genome if needed
if not args.bed_genome == args.extra_insight_genome:
    chosen_chain = chain_dict[args.bed_genome][args.extra_insight_genome]
    # Liftover bed file to EI co-ordinates and store in ei_1_bed
    tmp, ei_1_bed = sitewise_liftover(bed_local, chosen_chain, bed_dir,
                                          args.bed_genome, args.extra_insight_genome)
    sites_after_lo = count_sites(ei_1_bed)
    print(f"Sites after EI liftover({args.bed_genome} -> {args.extra_insight_genome}): {sites_after_lo}")
    # Successfuly mapped regions in source co-ordinates not needed
    os.remove(tmp)
else:
    ei_1_bed = bed_local

# Filter out sites not in the extraINSIGHT database
print("Filtering out sites not in the extraINSIGHT database")
ei_2_bed = os.path.join(bed_dir, f"ei_filtered.{args.extra_insight_genome}.bed")
cmd_ei_filter = f"/usr/bin/env bash -c \"zcat {ei_mut_cover_path} | bedops -i - <(zcat {ei_1_bed}) > {ei_2_bed}\""
# This fancy bit requred to keep error from being thrown when gunzip is interrupted when on of the input files terminates
# before the other. It's really annoying but harmless. Anyway this stops the pointless error
subprocess.call(args = cmd_ei_filter, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True)
# bash_call(cmd_ei_filter, True)
sites_after_filter = count_sites(ei_2_bed)
print(f"Sites after EI filtering: {sites_after_filter}")
os.remove(ei_1_bed)


# Liftover the extraINSIGHT bed file to the insight genome if needed
if not args.extra_insight_genome == args.insight_genome:
    chosen_chain = chain_dict[args.extra_insight_genome][args.insight_genome]
    ei_final_bed, i_bed = sitewise_liftover(ei_2_bed, chosen_chain, bed_dir,
                                            args.extra_insight_genome, args.insight_genome)
    sites_after_lo2 = count_sites(i_bed)
    print(f"Sites after I liftover({args.extra_insight_genome} -> {args.insight_genome}): {sites_after_lo2}")
    sites_after_lo2 = count_sites(ei_final_bed)
else:
    ei_final_bed = ei_2_bed
    i_bed = ei_2_bed

# Rename output (actually copy for robustness but whatever)
ei_path_final = os.path.join(bed_dir, "ei_" + args.extra_insight_genome + ".bed")
i_path_final = shutil.copyfile(i_bed, os.path.join(bed_dir, "i_" + args.insight_genome + ".bed"))
shutil.copyfile(ei_2_bed, ei_path_final)
shutil.copyfile(i_bed, i_path_final)

# Remove intermediate files
for i in {i_bed, ei_2_bed, ei_final_bed}:
    os.remove(i)

i_bed = os.path.abspath(i_path_final)
ei_bed = os.path.abspath(ei_path_final)

################### Done filtering bed files! ######################

# Setup insight analysis directories
extra_dir = os.path.abspath(os.path.join(args.out_dir, "ExtRaINSIGHT"))
make_directory(extra_dir)
insight_dir = os.path.abspath(os.path.join(args.out_dir, "INSIGHT"))
make_directory(insight_dir)

cwd = os.getcwd()
# Run ExtRaINSIGHT
print("Running ExtRaINSIGHT...")
os.chdir('../../extraINSIGHT/')
cmd_extra = f"""./ExtRaINSIGHT.py -b {ei_path_final} -o {extra_dir} \
-r {ei_mut_rate_path}"""
bash_call(cmd_extra)

os.chdir(cwd)

# Run INSIGHT
print("Running INSIGHT2")
cmd_insight = f"../../INSIGHT2/INSIGHT2.py -b {i_path_final} -o {insight_dir}"
bash_call(cmd_insight)
