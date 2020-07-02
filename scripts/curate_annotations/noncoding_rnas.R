library(GenomeInfoDb)
library(GenomicRanges)
library(data.table)
library(stringr)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

## Gets bases 2-7
seed_region <- function(gr) {
    start(gr[strand(gr) == "+"]) <- start(gr[strand(gr) == "+"]) + 1
    end(gr[strand(gr) == "+"]) <- pmin(start(gr[strand(gr) == "+"]) + 5,
                                       end(gr[strand(gr) == "+"]))
    end(gr[strand(gr) == "-"]) <- end(gr[strand(gr) == "-"]) - 1
    start(gr[strand(gr) == "-"]) <- pmax(start(gr[strand(gr) == "-"]),
                                       end(gr[strand(gr) == "-"]) - 5)
    return(gr)    
}

############################
## Micro-RNA seed regions
############################

## Out directories
out_data <- "../../data/grch37/targetscan_v2.7"
out_anno <- "../../data/grch37/annotation_bed/mirna"
dir.create(out_data, FALSE, TRUE)
dir.create(out_anno, FALSE, TRUE)

## Download url from targetscan v7.2 on hg19
ts_url <- "http://www.targetscan.org/vert_72/vert_72_data_download/All_Target_Locations.hg19.bed.zip"

## File locations
archive_file <-file.path(out_data, basename(ts_url))
download.file(ts_url, archive_file)
unzip(archive_file, exdir = out_data)
unlink(archive_file)

for(d in dir(out_data, "Targets_CS_pctiles.hg19.*.bed$", full.names = T)) {
    R.utils::gzip(d, overwrite = T)
}

## Note: Conservation status is relative to mammals
header <- c("chrom","start","end","name","score","strand")
data_list = list()
for(d in dir(out_data, "Targets_CS_pctiles.hg19.*.bed.gz$", full.names = TRUE)) {
    data <- fread(d)[,1:6]
    setnames(data, colnames(data), header)
    file_info <- str_match(d, ".+\\.(.+)\\.(.+)\\.bed.gz")
    data[,mir_family_conservation := file_info[1,2]]
    data[,site_conservation := file_info[1,3]]
    data_list[[d]] <- data
}

## Export bed files
for(b in seq_along(data_list)) {
    out_name <- paste0("targetscan_2.7_",
                       data_list[[b]]$mir_family_conservation[1], "_",
                       data_list[[b]]$site_conservation[1], "_seed_hg19.bed.gz")
    gr <- with(data_list[[b]], GRanges(chrom, IRanges(start, end - 1, names = name),
                                       strand = strand))
    gr <- gr[width(gr) <= 8]
    gr_seed <- reduce(seed_region(gr))
    export_bed(gr_seed, out_name, out_anno)
}


##################################################
## split 3' UTR into micro RNA and non-micro RNA
##################################################
## Out directories
out_data <- "../../data/grch37/targetscan_v2.7"
out_anno <- "../../data/grch37/annotation_bed/three_utr_decomposition"
dir.create(out_anno, FALSE, TRUE)

## Setup filters
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)

## Load txdb
txdb_37 <- load_txdb("GRCh37")

## Extract a variety of regions
cds_gr <- unlist(cdsBy(txdb_37, by = "tx", filter = filters))
three_utr_gr <- unlist(threeUTRsByTranscript(txdb_37, filter = filters))
three_utr_gr <- setdiff(three_utr_gr, cds_gr, ignore.strand = TRUE)
GenomeInfoDb::seqlevelsStyle(three_utr_gr) <- "ensembl"
three_utr_gr <- sort(sortSeqlevels(three_utr_gr))

## Note: Conservation status is relative to mammals
header <- c("chrom","start","end","name","score","strand")
data_list = list()
for(d in dir(out_data, "Targets_CS_pctiles.hg19.*.bed.gz$", full.names = TRUE)) {
    data <- fread(d)[,1:6]
    setnames(data, colnames(data), header)
    file_info <- str_match(d, ".+\\.(.+)\\.(.+)\\.bed.gz")
    data[,mir_family_conservation := file_info[1,2]]
    data[,site_conservation := file_info[1,3]]
    data_list[[d]] <- data
}

## Merge all microRNAs into one list and remove sex chromosomes
all_mirna <- rbindlist(data_list)
gr <- with(all_mirna, GRanges(chrom, IRanges(start, end - 1, names = name),
                              strand = strand))
GenomeInfoDb::seqlevelsStyle(gr) <- "ensembl"
gr <- keepSeqlevels(gr, as.character(1:22), pruning.mode="coarse")
gr <- sort(sortSeqlevels(gr))

## Extract seed regions and reduce
gr_seed <- seed_region(gr)
gr_seed <- sort(reduce(gr_seed))

## Partition utr into seed and non-seed
no_seed_utr <- setdiff(three_utr_gr, gr_seed)
seed_utr <- intersect(gr_seed, three_utr_gr)

export_bed(no_seed_utr, file_name = "utr_no_seed.hg19.bed.gz", path = out_anno)
export_bed(seed_utr, file_name = "utr_seed.hg19.bed.gz", path = out_anno)

#################
## Noncoding RNAs
#################

## Out directories
out_anno <- "../../data/grch38/annotation_bed/noncoding_rna"
dir.create(out_anno, FALSE, TRUE)

## Setup filters
seqnames_filter <- SeqNameFilter(as.character(1:22))
cds_gr <- reduce(unlist(cdsBy(txdb_38, by = "tx", filter = seqnames_filter)))

## Load txdb
txdb_38 <- load_txdb("GRCh38")
tx_biotypes <- c("lncRNA", "miRNA", "snoRNA", "snRNA")

for(txb in tx_biotypes) {
    out_name <- paste0(txb, ".bed.gz")
    txbiotype_filter <- TxBiotypeFilter(txb)
    filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)
    tx <- transcripts(txdb_38, filter = filters)
    tx <- sort(sortSeqlevels(setdiff(tx, cds_gr)))
    export_bed(tx, file_name = out_name, path = out_anno)
}
