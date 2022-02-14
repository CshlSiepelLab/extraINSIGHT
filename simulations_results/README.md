This directory contains scripts and parameters for generating simulations and results found in:  

Extreme purifying selection against point mutations in the human genome
Noah Dukler, Mehreen R. Mughal, Ritika Ramani, Yi-Fei Huang, Adam Siepel
bioRxiv 2021.08.23.457339; doi: https://doi.org/10.1101/2021.08.23.457339

The script binom.py is used to generate the results found in Figures S1 and S2. An example of the neutral mutation model file is included in the directory and can be run as follows:

python binom.py neutral_model_example 

This will generate a count of rare variants as predicted by our Bernoulli model as well as the observed rare variants for each line in the neutral mutation model file.

The directory "DFEs" contains simulation parameters for results shown in Figure S9 and Table S2. Specifically, each line contains a selection coefficient and number of sites used to generate the "Simulated DFE". The simulator used to generate these results in described in [1] and was shared with us by the authors. There are four files in this directory, each named according to the genomic element we are aiming to simulate.

[1] Donate Weghorn, Daniel J Balick, Christopher Cassa, Jack A Kosmicki, Mark J Daly, David R Beier, Shamil R Sunyaev
    Applicability of the Mutation–Selection Balance Model to Population Genetics of Heterozygous Protein-Truncating Variants in Humans
    Molecular Biology and Evolution, Volume 36, Issue 8, August 2019, Pages 1701–1710, https://doi.org/10.1093/molbev/msz092
