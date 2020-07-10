#!/usr/bin/env bash

## A script designed to take a chain file and create a bed file showing its coverage on the target genome
## where chain has query --> target directionality. Places output file in same directory as source file

bname=`basename $1`
target=$2
query=$3
oname=`readlink -e $1`
oname=`dirname ${oname}`/${target}_covered_by_${query}.bed.gz

if [[ ${bname} == *.gz ]]; then
    zcat $1 | awk 'length {if($1 ~ "chain"){chrom=$8;qChrLen=$9;qStrand=$10; if(qStrand=="+"){start=$11} else{start=qChrLen-$11}} else{if(qStrand=="+"){print chrom,start,start+$1;start=start+$1+$3}; if(qStrand=="-"){print chrom,start-$1,start;start=start-$1-$3}}}' | sort-bed - | bedops -m - | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | gzip -c > ${oname}
else    
    awk 'length {if($1 ~ "chain"){chrom=$8;qChrLen=$9;qStrand=$10; if(qStrand=="+"){start=$11} else{start=qChrLen-$11}} else{if(qStrand=="+"){print chrom,start,start+$1;start=start+$1+$3}; if(qStrand=="-"){print chrom,start-$1,start;start=start-$1-$3}}}' $1 | sort-bed - | bedops -m - | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | gzip -c > ${oname}
fi
