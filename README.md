# extraINSIGHT

## Description
A model for estimating the fraction of sites under strong negative selection by observing a depletion of rare variants within a set of genomic regions with respect to a neutral mutation model. This repository contains the pipeline for extacting the relevant data and fitting the neutral mutation model as well as the core extraINSIGHT model itself. This pipeline is available as a webserver from the Siepel Lab and can be found at http://compgen.cshl.edu/extrainsight

The accompanying manuscript is available at:

Extreme purifying selection against point mutations in the human genome Noah Dukler, Mehreen R. Mughal, Ritika Ramani, Yi-Fei Huang, Adam Siepel bioRxiv 2021.08.23.457339; doi: https://doi.org/10.1101/2021.08.23.457339

## Dependencies
- bedops (>= 2.4.2)
- bedtools (>= 2.22.0)
- conda (>= 4.7.10)
- R (>= 3.5.1)
- perl (>= 5.16.3)
