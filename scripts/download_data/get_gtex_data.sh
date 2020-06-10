#! /usr/bin/env bash

outdir=../../data/grch38/gtex
mkdir -p ${outdir}

wget -b -O ${outdir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz


