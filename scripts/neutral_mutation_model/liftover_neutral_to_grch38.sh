####################
## This script lifts over the neutral regions from hg19 to hg38.
## This is not used by the snakemake pipeline it must be run seperately.
####################

proj_root=/sonas-hs/siepel/nlsas/data/projects/extraINSIGHT
chainfile=${proj_root}/data/grch38/liftover/hg19ToHg38.over.chain.gz
neutral_dir=${proj_root}/data/grch38/neutral_regions

mkdir -p ${neutral_dir}

liftOver ${proj_root}/data/grch37/neutral_regions/hg19_neutral_region.bed ${chainfile} ${neutral_dir}/grch38_neutral_region.bed ${neutral_dir}/grch38_neutral_region.unmapped.bed
