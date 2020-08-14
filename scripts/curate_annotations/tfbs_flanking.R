library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/tfbs"
dir.create(out_path, FALSE, TRUE)

txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)

## Extract a variety of regions
promoters_gr <- reduce(promoters(txdb, filter = filters, upstream = 2000, downstream = 200))

## URL to ensembl regulatory build
## "ftp://ftp.ensembl.org/pub/release-100/regulation/homo_sapiens/MotifFeatures/Homo_sapiens.GRCh38.motif_features.gff.gz"
reg_gff <- "~/Projects/extraINSIGHT/data/grch38/ensembl_regulatory_build_100/Homo_sapiens.GRCh38.motif_features.gff.gz"
tmp <- tempdir()
export_bed(promoters_gr, file_name = "promoters.bed.gz", path = tmp, style = "ensembl")
query_bed <- file.path(tmp, "promoters.bed.gz")
out <- file.path(tmp, "out.gff")

## Query gff for tfbs in promoters
cmd <- paste("bash -c \' tabix", reg_gff, "-R <( zcat", query_bed, " | sort-bed - ) >", out, "\'")
system(cmd)
gr <- import(out, format = "GFF")
unlink(out)

## Switch to 0-based inclusive
start(gr) <- start(gr) - 1
end(gr) <- end(gr) - 1

## Get flanks
gr_flank <- reduce(flank(gr, width = 10), ignore.strand = T)

## Reduce gr for tfbs
gr <- reduce(gr, ignore.strand = T)

## Get disjoint sets of tfbs and flanking regions
gr_flank_final <- reduce(sort(setdiff(gr_flank, gr)))
gr_final <- reduce(sort(setdiff(gr, gr_flank)))

export_bed(gr_final, "tfbs.bed.gz", path = out_path)
export_bed(gr_flank_final, "tfbs_flank.bed.gz", path = out_path)
