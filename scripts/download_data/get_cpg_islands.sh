## This script downloads and processes CpG island annotations

# YAML parser from stackoverflow
# https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Read yaml variables into the local environment
grch38_yaml=../config_files/grch38_neutral_config.yml
eval $(parse_yaml $grch38_yaml "GRCH38_")

grch37_yaml=../config_files/grch37_neutral_config.yml
eval $(parse_yaml $grch37_yaml "GRCH37_")


grch38_outdir=${GRCH38_rootdir}/`dirname ${GRCH38_covariates_GC_COVARIATE_FILE}`
grch37_outdir=${GRCH37_rootdir}/`dirname ${GRCH37_covariates_GC_COVARIATE_FILE}`

# From https://www.biostars.org/p/236141/
# Derived from the table schema for this file, the first four columns are the island's
# genomic interval and name. The remaining columns are island length, number of CpGs in
# the island, the number of C and G in the island, the percentage of island that is CpG,
# the percentage of island that is C or G, and the ratio of observed(cpgNum) to
# expected(numC*numG/length) CpG in island. NOTE: 0-based start co-ordinates half open

# Get and process CPG islands for hg38
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, substr($0, index($0, $7)); }' \
   | sort-bed - \
   | gzip -c \
   > ${grch38_outdir}/cpgIslandExt.hg38.bed.gz

# and grch37
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, substr($0, index($0, $7)); }' \
   | sort-bed - \
   | gzip -c \
   > ${grch37_outdir}/cpgIslandExt.hg37.bed.gz
