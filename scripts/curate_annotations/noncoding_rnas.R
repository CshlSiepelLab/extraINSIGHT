library(GenomeInfoDb)
library(GenomicRanges)
library(data.table)
library(stringr)
# library(openxlsx)
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

## Extract transcripts
txdb <- load_txdb("GRCh37")
tx_gr <- transcripts(txdb)

## Remove all nonstandard and seq chromosomes
tx_gr <- keepSeqlevels(tx_gr, as.character(1:22), pruning.mode="coarse")

## Get only protein coding transcripts
filter_tx <- unlist(tx_gr[tx_gr$tx_biotype == "protein_coding",]$tx_id)
tx_gr <- tx_gr[tx_gr$tx_id%in% filter_tx]

## Extract a variety of regions
cds_gr <- unlist(cdsBy(txdb, by = "tx")[filter_tx])
three_utr_gr <- unlist(threeUTRsByTranscript(txdb))
three_utr_gr <- three_utr_gr[names(three_utr_gr) %in% filter_tx]
three_utr_gr <- setdiff(three_utr_gr, cds_gr, ignore.strand = TRUE)

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
gr <- sort(gr)

## Extract seed regions and reduce
gr_seed <- seed_region(gr)
gr_seed <- reduce(gr_seed)

## Partition utr into seed and non-seed
no_seed_utr <- setdiff(three_utr_gr, gr_seed)
seed_utr <- intersect(gr_seed, three_utr_gr)

export_bed(no_seed_utr, file_name = "utr_no_seed.hg19.bed.gz", path = out_anno)
export_bed(seed_utr, file_name = "utr_no_seed.hg19.bed.gz", path = out_anno)

#################
## lncRNAs
#################
## out_lnc_data <- "../../data/grch38/lncRNA"
## dir.create(out_lnc_data, F, T)
## ## Retrive lncRNA coordinates from RNACentral v15
## ## https://rnacentral.org/
## lnc_rna_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/15.0/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz"
## lnc_rna_bed <- fread(lnc_rna_url)

## ## Split into different types of ncRNAs
## table(lnc_rna_bed)


## ## Disease associated LNCrna from the lncRNA disease database
## ## http://www.rnanut.net/lncrnadisease/index.php/home
## lnc_disease_url <- "http://www.rnanut.net/lncrnadisease/static/download/all%20ncRNA-disease%20information.xlsx"
##                                         # lnc_rna_disease <- fread(lnc_disease_url)
## disease_xlsx_path <- file.path(out_lnc_data, "all_ncRNA-disease_information.xlsx")
## download.file(url = lnc_disease_url, destfile = disease_xlsx_path, mode="wb")
## disease_lnc_rna <- as.data.table(read.xlsx(disease_xlsx_path))
## disease_lnc_rna[Species == "Homo sapiens"][1,]
