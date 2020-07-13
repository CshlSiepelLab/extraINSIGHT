library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
source("load_txdb.R")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/phastcons_elements"
dir.create(out_path, FALSE, TRUE)

## Download from ucsc
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/phastConsElements100way.txt.gz"
bed <- fread(url)
bed_gr <- with(bed, GRanges(V2, IRanges(V3, V4 - 1)))

## Convert to ensembl style and remove non-autosomes
seqlevelsStyle(bed_gr) <- "Ensembl"
bed_gr <- keepSeqlevels(bed_gr, 1:22, pruning.mode = "coarse")

## By definition there are no CDS in phastCons elements so no need to remove them
## (I already double checked this anyway) ... so just export them
export_bed(reduce(bed_gr), "phastcons_elements_100way.bed.gz", path = out_path)


