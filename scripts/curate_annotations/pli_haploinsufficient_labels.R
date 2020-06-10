library(GenomicFeatures)
library(GenomicRanges)
library(data.table)

export_bed <- function(gr, file_name, path = ".", adjust_start = 0, adjust_end = 0, style = "UCSC",
                       sort_bed= TRUE) {
    dir.create(path, FALSE, TRUE)
    start(gr) <- start(gr) + adjust_start
    end(gr) <- end(gr) + adjust_end
    out_file <- file.path(path, file_name)
    GenomeInfoDb::seqlevelsStyle(gr) <- style
    if (sort_bed) {
        seqlevels(gr) <- sort(as.character(seqlevels(gr)))
        gr <- sort(gr)
    }
    if(tools::file_ext(file_name) == "gz") {
        con <- gzfile(out_file)
    } else {
        con <- gzfile(out_file)
    }
    rtracklayer::export.bed(gr, con)
}

## Ensembl release 99 = Gencode v?? but on grch37
## ensembl_gff <- "../../data/grch37/gtf/Homo_sapiens.GRCh37.87.gtf.gz"
## txdb <- makeTxDbFromGFF(ensembl_gff, format = "gtf")
## saveDb(txdb, file="grch37.sqlite")

txdb <- loadDb("grch37.sqlite")

## Load in transcript database
out_path <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/pli_annotations"
dir.create(out_path, FALSE, TRUE)

## Read in pLI file
pli_file <- "../../data/grch37/pli_v2.1.1/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
pli <- fread(cmd = paste("bgzip -d -c", pli_file))
setkey(pli, "chromosome")

## Label genes as above or below pLI threshold, keep only chromsomes 1-22
pli_out <-  pli[.(as.character(1:22)), .(transcript, pLI_09 = pLI >= 0.9)]
pli_hi <- pli_out[pLI_09 == TRUE,]$transcript
pli_lo <- pli_out[pLI_09 == FALSE,]$transcript

## Get exons that are solely in one class or the other
exons_hi <- exons(txdb, filter = list("TXNAME" = pli_hi))
exons_lo <- exons(txdb, filter = list("TXNAME" = pli_lo))
exons_shared <- intersect(exons_hi, exons_lo)
exons_hi <- setdiff(exons_hi, exons_shared)
exons_lo <- setdiff(exons_lo, exons_shared)

export_bed(exons_hi, file_name = "pli_high_exons.bed.gz", path = out_path)
export_bed(exons_lo, file_name = "pli_low_exons.bed.gz", path = out_path)

## Get cds that are solely in one class or the other
cds_hi <- cds(txdb, filter = list("TXNAME" = pli_hi))
cds_lo <- cds(txdb, filter = list("TXNAME" = pli_lo))
cds_shared <- intersect(cds_hi, cds_lo)
cds_hi <- setdiff(cds_hi, cds_shared)
cds_lo <- setdiff(cds_lo, cds_shared)

export_bed(cds_hi, file_name = "pli_high_cds.bed.gz", path = out_path)
export_bed(cds_lo, file_name = "pli_low_cds.bed.gz", path = out_path)

## Now 5'-UTRs
cds_gr <- cds(txdb)
five_utr_gr <- unlist(fiveUTRsByTranscript(txdb, use.names = TRUE))
five_utr_hi <- five_utr_gr[names(five_utr_gr) %in% pli_hi]
five_utr_lo <- five_utr_gr[names(five_utr_gr) %in% pli_lo]
five_utr_hi <- setdiff(five_utr_hi, cds_gr, ignore.strand = TRUE)
five_utr_lo <- setdiff(five_utr_lo, cds_gr, ignore.strand = TRUE)
five_utr_shared <- intersect(five_utr_hi, five_utr_lo)
five_utr_hi <- setdiff(five_utr_hi, five_utr_shared)
five_utr_lo <- setdiff(five_utr_lo, five_utr_shared)

export_bed(five_utr_hi, file_name = "pli_high_five_utr.bed.gz", path = out_path)
export_bed(five_utr_lo, file_name = "pli_low_five_utr.bed.gz", path = out_path)

## Now 3'-UTRs
cds_gr <- cds(txdb)
three_utr_gr <- unlist(threeUTRsByTranscript(txdb, use.names = TRUE))
three_utr_hi <- three_utr_gr[names(three_utr_gr) %in% pli_hi]
three_utr_lo <- three_utr_gr[names(three_utr_gr) %in% pli_lo]
three_utr_hi <- setdiff(three_utr_hi, cds_gr, ignore.strand = TRUE)
three_utr_lo <- setdiff(three_utr_lo, cds_gr, ignore.strand = TRUE)
three_utr_shared <- intersect(three_utr_hi, three_utr_lo)
three_utr_hi <- setdiff(three_utr_hi, three_utr_shared)
three_utr_lo <- setdiff(three_utr_lo, three_utr_shared)

export_bed(three_utr_hi, file_name = "pli_high_three_utr.bed.gz", path = out_path)
export_bed(three_utr_lo, file_name = "pli_low_three_utr.bed.gz", path = out_path)
