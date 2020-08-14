library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
library(data.table)
source("annotation_lib.R")
source("load_txdb.R")

txdb_38 <- load_txdb("GRCh38")
txdb_37 <- load_txdb("GRCh37")

#####################
## Get coding annotations to calculate genomewide averages for
## constraint on coding regions for both grch37 and grch38
#####################
out_path_38 <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/coding_genomewide"
out_path_37 <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/coding_genomewide"
dir.create(out_path_38, FALSE, TRUE)
dir.create(out_path_37, FALSE, TRUE)

seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)

cds_gr_37 <- unlist(cdsBy(txdb_37, filter = filters))
cds_gr_38 <- unlist(cdsBy(txdb_38, filter = filters))

strand(cds_gr_37) <- "*"
strand(cds_gr_38) <- "*"
cds_gr_37 <- reduce(cds_gr_37, ignore.strand = T)
cds_gr_38 <- reduce(cds_gr_38, ignore.strand = T)

export_bed(cds_gr_37, "coding.bed.gz", path = out_path_37)
export_bed(cds_gr_38, "coding.bed.gz", path = out_path_38)


#####################
## Get non-coding annotations to calculate genomewide averages for
## constraint on non-coding regions for both grch37 and grch38
#####################
out_path_38 <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/noncoding_genomewide"
out_path_37 <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/noncoding_genomewide"
dir.create(out_path_38, FALSE, TRUE)
dir.create(out_path_37, FALSE, TRUE)

noncode_37 <- gaps(cds_gr_37)
noncode_38 <- gaps(cds_gr_38)
noncode_37 <- noncode_37[strand(noncode_37) == "*"]
noncode_38 <- noncode_38[strand(noncode_38) == "*"]

export_bed(sort(reduce(noncode_37)), "noncoding.bed.gz", path = out_path_37)
export_bed(sort(reduce(noncode_38)), "noncoding.bed.gz", path = out_path_38)
