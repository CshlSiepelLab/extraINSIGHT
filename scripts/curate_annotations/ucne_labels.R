library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

## UNCE annotations from UNCEbase
## https://academic.oup.com/nar/article/41/D1/D101/1057253
out_dir <- "../../data/grch37/annotation_bed/ucne"
dir.create(out_dir, FALSE, TRUE)
ucne <- fread("ftp://ccg.vital-it.ch/UCNEbase/ucnes/hg19_UCNE_coord.bed")

txdb <- load_txdb("GRCh37")
seqnames_filter <-  AnnotationFilterList(SeqNameFilter(as.character(1:22)))
filters <- c(seqnames_filter)

## Write ucne to bed file
fwrite(ucne[,.(V1,V2,V3)], file = file.path(out_dir, "ucne_all.bed.gz"), sep = "\t", col.names = F)

## Convert UCNE to GRanges
ucne_gr <- with(ucne, GRanges(V1, IRanges(V2, V3 - 1)))

## Get intronic regions from database
five_prime_utr <- reduce(unlist(fiveUTRsByTranscript(txdb, filter = filters)))
GenomeInfoDb::seqlevelsStyle(five_prime_utr) <- "UCSC"
introns_gr<- reduce(unlist(intronsByTranscript(txdb, filter = filters)))
GenomeInfoDb::seqlevelsStyle(introns_gr) <- "UCSC"
transcripts_gr <- reduce(transcripts(txdb, filter = filters))
GenomeInfoDb::seqlevelsStyle(transcripts_gr) <- "UCSC"

## Get UCNE that overlap with introns, then remove regions that are in or near 5'UTR
## because of known mutation rate issues
ucne_intra_tx <- intersect(ucne_gr, introns_gr, ignore.strand = TRUE)
                         
## Get UCNE are not in genes or their promoter regions
ucne_inter_tx <- setdiff(ucne_gr, reduce(c(transcripts_gr,
                                    promoters(transcripts_gr,
                                              upstream = 2000,
                                              downstream = 200)
                                    )), ignore.strand = TRUE)

export_bed(reduce(ucne_intra_tx), "ucne_intra_transcript.bed.gz", out_dir)
export_bed(reduce(ucne_inter_tx), "ucne_inter_transcript.bed.gz", out_dir)



