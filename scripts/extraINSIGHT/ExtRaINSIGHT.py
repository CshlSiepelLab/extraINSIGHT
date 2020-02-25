#! /usr/bin/env python

from __future__ import division, print_function
import sys
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import chi2
import argparse
import ExtRaINSIGHT_fun

#
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="input_file", nargs=1, type=str,
                    required=True, help="input file")
args = parser.parse_args()

obs, neutral_rate = read_data(args.input_file[0])

if obs is None:
    print("Error: Unable to estimate selection")
    sys.exit()
    
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
