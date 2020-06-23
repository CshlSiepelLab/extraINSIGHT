#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import chi2
import argparse
import os
import pandas as pd
import gzip
import io


# Specify parser
parser = argparse.ArgumentParser()
parser.add_argument("-b","--bed", dest="bed", nargs=1, type=str,
                    required=True, help="bed file containing query regions")
parser.add_argument("-r","--mutation-rate", dest="mutation_rate", nargs=1, type=str,
                    required=True, help="bed file containing per allele mutation rates")
parser.add_argument("-o","--out-dir", dest="out_dir", nargs=1, type=str,
                    required=True, help="results output directory")

args = parser.parse_args()

args.out_dir = args.out_dir[0]
args.bed = args.bed[0]

# File paths
anno_bed_path = os.path.join(args.out_dir, "annotated.bed.gz")
final_bed_path = os.path.join(args.out_dir, "final.bed.gz")
ei_out_path = os.path.join(args.out_dir, "strong_selection_estimate.txt")

cmd_anno = f"""tabix {args.mutation_rate[0]} -R {args.bed[0]} | intersectBed -a - -b {args.bed[0]} -sorted |\
gzip -c > {anno_bed_path}"""
os.system(cmd_anno)
cmd_final_bed=f"zcat {anno_bed_path} | cut -f1-3 | bedops -m - | gzip -c > {final_bed_path}"
os.system(cmd_final_bed)

df = pd.read_csv(anno_bed_path, sep = "\t", header = None, usecols = [4, 6])

# # Seperate into two seperate vectors
obs = df.iloc[:, 0].to_numpy()
neutral_rate = df.iloc[:, 1].to_numpy()

# If no sites that met filters in region, kill job
if obs is None:
    raise ValueError("Error: Unable to estimate selection due to lack of sites that pass quality thresholds")
    
# Computes the negative likelihood of the observed number of mutations under the
# neutral rates
def mutation_negative_log_likelihood(rate, obs, neutral_rate):
    assert obs.size == neutral_rate.size, "Unmatched number of sites"
    expected_rate = rate * neutral_rate
    probs = np.where(obs , expected_rate, 1. - expected_rate)
    nll = -np.sum(np.log(probs))
    return(nll)

print("Fitting ExtRaINSIGHT model ...") 
# Set heuristic upper bound for optimizer
max_rate = np.min(1. / neutral_rate) - 1.e-3

# fit alternative model
res = minimize(mutation_negative_log_likelihood,
               1., args=(obs, neutral_rate),
               method="L-BFGS-B", bounds=[(1.e-3, max_rate)])
relative_rate = res.x[0]
nll_1 = res.fun

# fit null model
nll_0 = mutation_negative_log_likelihood(1., obs, neutral_rate)
    
# likelihood ratio test
statistic = 2. * (nll_0 - nll_1)
# log_p_value = chi2.logsf(statistic, 1) / np.log(10.)

# the likelihood ratio statistic follows a chi-square distribution
p_value = chi2.sf(statistic, 1)
    
# calculate curvature
curvature = sp.misc.derivative(mutation_negative_log_likelihood,
                               relative_rate, dx=1.e-6,
                               args=(obs, neutral_rate),
                               n=2)

# print output
nm = ['strong_selection', 'p', 'lr_stat', 'num_possible_mutations', 'curvature', 'se']
val = [1. - relative_rate, p_value, statistic, obs.size, curvature, np.sqrt(1. / curvature)]
df_out = pd.DataFrame(list(zip(nm, val)), columns =['name', 'val']) 

df_out.to_csv(ei_out_path, index=False, sep='\t')
