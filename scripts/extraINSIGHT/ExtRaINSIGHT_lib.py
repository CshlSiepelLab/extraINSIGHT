# Contains functions used by ExtRaINSIGHT
import numpy as np

# Computes the negative likelihood of the observed number of mutations under the
# neutral rates
def mutation_negative_log_likelihood(rate, obs, neutral_rate):
    assert obs.size == neutral_rate.size, "Unmatched number of sites"
    expected_rate = rate * neutral_rate
    probs = np.where(obs == 1., expected_rate, 1. - expected_rate)
    nll = -np.sum(np.log(probs))
    # print("x = {}; y = {}".format(sel, nll))
    return(nll)

# Reads in bed file
def read_data(in_file):
    # load data to array
    data = np.loadtxt(in_file, delimiter="\t", usecols=(4, 6))

    obs = data[:, 0]
    neutral_rate = data[:, 1]

    return obs, neutral_rate
