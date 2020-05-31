library(data.table)
library(stringr)
library(stringdist)
library(ggplot2)

## Read in and parse context table
freq_table=fread("~/Projects/extraINSIGHT/results/grch37/gnomad_v2.1/mutation_model/mutation_frequency_table.txt")
freq_table[,`:=`(context=str_split(freq_table$V1,"-",simplify=TRUE)[,1],
                 alt=str_split(freq_table$V1,"-",simplify=TRUE)[,2])]
freq_table[,ref:=str_sub(context,4,4)]
freq_table[,reduced_context:=paste0(str_sub(context,1,3),str_sub(context,5,7))]
freq_table[,mutation_label:=paste0(ref,"->",alt)]

## Compute distances between contexts
context_distances = with(freq_table,stringdistmatrix(unique(reduced_context),unique(reduced_context)),
    method="lv",nthread =4)
colnames(context_distances) = unique(freq_table$reduced_context)
rownames(context_distances) = unique(freq_table$reduced_context)
## Produce heirarchical clustering of string contexts
temp = dist(context_distances)
context_ord = rownames(context_distances)[hclust(d=context_distances)$ord]
