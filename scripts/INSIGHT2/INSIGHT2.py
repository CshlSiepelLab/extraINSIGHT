#! /usr/bin/env python

import argparse
import os
import shlex
import shutil

# Path to the INSIGHT package
insight2_db = "/local/storage/data/extraINSIGHT/no-backup/data/INSIGHT2/Insight2DB"
insight2 = "/local/projects/extraINSIGHT/bin/Insight2"
# insight2 = "/sonas-hs/siepel/nlsas/data/projects/extraINSIGHT/bin/Insight2"

# Set textwap width for parser
os.environ['COLUMNS'] = "90"

# Specify parser
parser = argparse.ArgumentParser()
parser.add_argument("-b","--bed", dest="bed", type=str,
                    required=True, help="bed file containing query regions")
parser.add_argument("-o","--out-dir", dest="out_dir", type=str,
                    default = ".", help="results output directory (default = '.')")
args = parser.parse_args()

print(args)

# Create output directory
if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)
    print("Directory " , args.out_dir ,  " Created ")
else:    
    print("Directory " , args.out_dir ,  " already exists")

# Expand out paths
args.bed = os.path.abspath(os.path.realpath(args.bed))
args.out_dir = os.path.abspath(os.path.realpath(args.out_dir))

# Get corename of bed file
data_id = os.path.splitext(os.path.basename(args.bed))[0]

# Copy the query bed file to the output directory
shutil.copy(args.bed, args.out_dir)

# remove bed suffix for insight2
bed_string = os.path.join(args.out_dir, data_id)

# Create command line call
cmd = f"{insight2} {insight2_db} -fin {bed_string} -qmap"
print(cmd)

# Run INSIGHT2, it will automatically create output in the same directory as the bed file
# which in this case is the output directory
os.system(cmd)
