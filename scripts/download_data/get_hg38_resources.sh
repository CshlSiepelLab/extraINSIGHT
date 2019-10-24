proj_root=/sonas-hs/siepel/nlsas/data/projects/extraINSIGHT
fa_dir=${proj_root}/data/grch38/fa

mkdir -p ${fa_dir}

# Download the genome sequence
wget -O ${fa_dir}/grch38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Generate the chromosome size file
faSize -detailed -tab ${fa_dir}/grch38.fa.gz > ${fa_dir}/grch38.chrom.sizes
