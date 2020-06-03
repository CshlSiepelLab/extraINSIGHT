library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

## Ensembl release 99 = Gencode v33
## txdb <- makeTxDbFromEnsembl(organism="Homo sapiens", release = 99)
## saveDb(txdb, file="grch38.sqlite")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/"

args = list()
args$txdb = "grch38.sqlite"
args$out_dir

txdb <- loadDb(args$txdb)

export_bed <- function(gr, file_name, path = ".", adjust_start = 0, adjust_end = 0, style = "UCSC") {
    dir.create(path, FALSE, TRUE)
    start(gr) <- start(gr) + adjust_start
    end(gr) <- end(gr) + adjust_end
    out_file <- file.path(path, file_name)
    GenomeInfoDb::seqlevelsStyle(gr) <- style
    if(tools::file_ext(file_name) == "gz") {
        con <- gzfile(out_file)
    } else {
        con <- gzfile(out_file)
    }
    rtracklayer::export.bed(gr, con)
}

## Extract transcripts
tx_gr <- transcripts(txdb, columns=c("GENEID", "TXTYPE", "TXNAME", "TXID"))

## Remove all nonstandard and seq chromosomes
tx_gr <- keepSeqlevels(tx_gr, as.character(1:22), pruning.mode="coarse")

## Get only protein coding transcripts
filter_tx <- unlist(tx_gr[tx_gr$TXTYPE == "protein_coding",]$TXID)
tx_gr <- tx_gr[tx_gr$TXID %in% filter_tx]

## Extract a variety of regions
exon_gr <- exons(txdb, filter = list("TXID" = filter_tx))
cds_gr <- cds(txdb, filter = list("TXID" = filter_tx))
five_utr_gr <- unlist(fiveUTRsByTranscript(txdb, use.names = FALSE))
five_utr_gr <- five_utr_gr[names(five_utr_gr) %in% filter_tx]
three_utr_gr <- unlist(threeUTRsByTranscript(txdb, use.names = FALSE))
three_utr_gr <- three_utr_gr[names(three_utr_gr) %in% filter_tx]
promoters_gr <- promoters(txdb, filter = list("TXID" = filter_tx), upstream = 2000, downstream = 200)


## Subtract 1 off from the start b/c ensembl is 1-based
export_bed(exon_gr, "exon_gencode33.bed.gz", path = out_path)
export_bed(cds_gr, "cds_gencode33.bed.gz", path = out_path)
export_bed(five_utr_gr, "five_utr_gencode33.bed.gz", path = out_path)
export_bed(three_utr_gr, "three_utr_gencode33.bed.gz", path = out_path)
export_bed(promoters_gr, "promoters_gencode33.bed.gz", path = out_path)

## Get introns and remove all exonic regions (ignoring strand)
introns_gr <- unlist(intronsByTranscript(txdb, use.names = FALSE))
introns_gr <- introns_gr[names(introns_gr) %in% filter_tx]
introns_gr <- setdiff(introns_gr, exon_gr, ignore.strand = TRUE)

export_bed(introns_gr, "introns_gencode33.bed.gz", path = out_path)

## Splice junctions
## for discussion, c.f.  [gaps does not work on a GRangesList]
## (https://github.com/Bioconductor/GenomicRanges/issues/17)
introns_tmp <- unlist(intronsByTranscript(txdb, use.names = FALSE))
introns_tmp <- introns_tmp[names(introns_tmp) %in% filter_tx]
splice_gr <- c(GRanges(seqnames(introns_tmp), IRanges(start(introns_tmp), start(introns_tmp) + 1)),
               GRanges(seqnames(introns_tmp), IRanges(end(introns_tmp) - 1, end(introns_tmp))))
splice_gr <- reduce(sort(splice_gr))

export_bed(splice_gr, "splice_gencode33.bed.gz", path = out_path)
