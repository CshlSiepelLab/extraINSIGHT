library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/gene_annotations"
dir.create(out_path, FALSE, TRUE)

txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)

## Extract transcripts
tx_gr <- transcripts(x = txdb, filter = filters)

## Extract a variety of regions
exon_gr <- exons(txdb, filter = filters)
cds_gr <- unlist(cdsBy(txdb, filter = filters))
five_utr_gr <- unlist(fiveUTRsByTranscript(txdb, filter = filters))
five_utr_gr <- setdiff(five_utr_gr, cds_gr, ignore.strand = TRUE)
three_utr_gr <- unlist(threeUTRsByTranscript(txdb, filter = filters))
three_utr_gr <- setdiff(three_utr_gr, cds_gr, ignore.strand = TRUE)
promoters_gr <- promoters(txdb, filter = filters, upstream = 2000, downstream = 200)

## Subtract 1 off from the start b/c ensembl is 1-based
export_bed(reduce(exon_gr), "exon_gencode33.bed.gz", path = out_path)
export_bed(reduce(cds_gr), "cds_gencode33.bed.gz", path = out_path)
export_bed(reduce(five_utr_gr), "five_utr_gencode33.bed.gz", path = out_path)
export_bed(sort(reduce(three_utr_gr)), "three_utr_gencode33.bed.gz", path = out_path)
export_bed(sort(reduce(promoters_gr)), "promoters_gencode33.bed.gz", path = out_path)

## Get introns and remove all exonic regions (ignoring strand)
all_exon_gr <- reduce(exons(txdb)) ## use all exons regardless of txtype to filter later
introns_gr <- unlist(intronsByTranscript(txdb, filter = filters))
introns_gr <- setdiff(introns_gr, all_exon_gr, ignore.strand = TRUE)

export_bed(sort(reduce(introns_gr)), "introns_gencode33.bed.gz", path = out_path)

## Splice junctions
## for discussion, c.f.  [gaps does not work on a GRangesList]
## (https://github.com/Bioconductor/GenomicRanges/issues/17)
introns_tmp <- unlist(intronsByTranscript(txdb, filter = filters))
splice_gr <- c(GRanges(seqnames(introns_tmp), IRanges(start(introns_tmp), start(introns_tmp) + 1)),
               GRanges(seqnames(introns_tmp), IRanges(end(introns_tmp) - 1, end(introns_tmp))))
splice_gr <- reduce(sort(splice_gr))

export_bed(reduce(splice_gr), "splice_gencode33.bed.gz", path = out_path)
