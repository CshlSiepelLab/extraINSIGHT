proj_root=/sonas-hs/siepel/nlsas/data/projects/extraINSIGHT
fa_dir=${proj_root}/data/grch38/fa
liftover_dir=${proj_root}/data/grch38/liftover

mkdir -p ${fa_dir}
mkdir -p ${liftover_dir}

# Download the genome sequence
wget -O ${fa_dir}/grch38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Generate the chromosome size file
faSize -detailed -tab ${fa_dir}/grch38.fa.gz > ${fa_dir}/grch38.chrom.sizes

# Download hg19 to hg38 chain to liftover neutral regions
wget -O ${liftover_dir}/hg19ToHg38.over.chain.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
