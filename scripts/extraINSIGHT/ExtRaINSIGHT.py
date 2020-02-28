#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import chi2
import argparse
import subprocess
import shlex


# Specify parser
parser = argparse.ArgumentParser()
parser.add_argument("-b","--bed", dest="bed", nargs=1, type=str,
                    required=True, help="bed file containing query regions")
parser.add_argument("-r","--mutation-rate", dest="mutation_rate", nargs=1, type=str,
                    required=True, help="bed file containing per allele mutation rates")
parser.add_argument("--max-sites", dest="max_sites", nargs=1, type=int,
                    required=False, help=
                    """sets the maximum of sites to be considered by ExtRaINSIGHT. If too many
                sites are specified they will be randomly subsampled to this value."""
)
parser.add_argument("-o","--out-dir", dest="out_dir", nargs=1, type=str,
                    required=True, help="results output directory")
parser.add_argument("--exclude-chr", dest="exclude_chr", nargs="+", type=str,
                    required=False, help="chromosomes to exclude from the analysis")

args = parser.parse_args()
# Make excluded chromosomes single string for passing
excl = "|".join(args.exclude_chr)
# Create snakemake call
snake_cmd = f"""snakemake --snakefile ExtRaINSIGHT.smk --config bed={args.bed[0]} max_sites={args.max_sites[0]} \
mutation_rate={args.mutation_rate[0]} out_dir={args.out_dir[0]} exclude_chr={excl}"""

## Call snakemake call and parse output to 2-D numpy array
split_snake = shlex.split(snake_cmd)
arr = np.loadtxt(subprocess.Popen(split_snake, stdout=subprocess.PIPE).communicate()[0].splitlines())

# Seperate into two seperate vectors
obs = arr[:,0]
neutral_rate = arr[:,1]

# If no sites that met filters in region, kill job
if obs is None:
    print("Error: Unable to estimate selection due to lack of sites that pass quality thresholds")
    sys.exit()

# Computes the negative likelihood of the observed number of mutations under the
# neutral rates
def mutation_negative_log_likelihood(rate, obs, neutral_rate):
    assert obs.size == neutral_rate.size, "Unmatched number of sites"
    expected_rate = rate * neutral_rate
    probs = np.where(obs == 1., expected_rate, 1. - expected_rate)
    nll = -np.sum(np.log(probs))
    # print("x = {}; y = {}".format(sel, nll))
    return(nll)

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
    
# print output
print(("strong_selection = {}\npvalue = {}\n"
       + "likelihood_ratio_statistic = {}\n"
       + "num_of_possible_mutations = {}")
      .format(1. - relative_rate, p_value, statistic,
              obs.size))

# calculate curvature
curvature = sp.misc.derivative(mutation_negative_log_likelihood,
                               relative_rate, dx=1.e-6,
                               args=(obs, neutral_rate),
                               n=2)
print("curvature = {}".format(curvature))
print("standard_error = {}".format(np.sqrt(1. / curvature)))
