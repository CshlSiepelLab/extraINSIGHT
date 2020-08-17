library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(GenomicFeatures)
source("annotation_lib.R")
source("load_txdb.R")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/micro_rnas"
dir.create(out_path, FALSE, TRUE)

## Load txdb for some filtering of flanking regions
txdb <- load_txdb("GRCh38")
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(txbiotype_filter)
exon_gr <- exons(txdb, filter = filters)
exon_gr <- reduce(exon_gr)
seqlevelsStyle(exon_gr) <- "ucsc"

## Download required files - mirbase version 22
gff_url <- "ftp://mirbase.org/pub/mirbase/22.1/genomes/hsa.gff3"
ts_url <- "http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip"
mirbase_gff <- import(gff_url)
temp <- tempfile()
download.file(ts_url,temp)
targetscan_df <- as.data.table(read.delim(unz(temp, "miR_Family_Info.txt")))
unlink(temp)

## Select human miRNAs
targetscan_df <- targetscan_df[targetscan_df$Species.ID == 9606, ]

## highly conserved (2), conserved (1), poorly conserved but confidently annotated (0)
## filter out poorly conserved and possibly misannotated as a miRNA (-1)
## http://www.targetscan.org/faqs.Release_7.html
targetscan_df <- targetscan_df[Family.Conservation. != -1 & MiRBase.Accession != ""]
mir_id_cons <- targetscan_df[, .(MiRBase.Accession, Family.Conservation.)]
setnames(mir_id_cons, colnames(mir_id_cons), c("ID", "conservation"))
highly_conserved <- mir_id_cons[mir_id_cons$conservation == 2]$ID
conserved <- mir_id_cons[mir_id_cons$conservation == 1]$ID
poorly_conserved <- mir_id_cons[mir_id_cons$conservation == 0]$ID

## separate pri-miRNA and miRNA matures
primary_mir <- mirbase_gff[mirbase_gff$type == "miRNA_primary_transcript"]
mature_mir <- mirbase_gff[mirbase_gff$type == "miRNA"]

## generate grouping for pri-miRNA and miRNA matures by the mir-rna that the mature transcript
## was derived from
primary_mir$group <- primary_mir$ID
mature_mir$group <- mature_mir$Derives_from

## Subset miRNAs by those present in targetscan for analysis
mature_mir <- mature_mir[mature_mir$ID %in% targetscan_df$MiRBase.Accession]
primary_mir <- primary_mir[primary_mir$group %in% mature_mir$group]

## separate miRNA matures into seed and non-seed
mcols(mature_mir) <- merge.data.frame(mcols(mature_mir), mir_id_cons,
                                      by = "ID", sort = FALSE, all.x = TRUE)
mature_mir$score <- mature_mir$conservation

## Get the first 8 bases of the range as the seed
mature_mir_seed <- promoters(mature_mir, 0, 8)
## Regions that do not overlap with that are not the seed region
mature_mir_nonseed <- setdiff(mature_mir, mature_mir_seed)
mcols(mature_mir_nonseed) <- mcols(mature_mir)[subjectHits(findOverlaps(mature_mir_nonseed, mature_mir)),]

## Get only conserved
mature_mir_seed_conserved <- mature_mir_seed[mature_mir_seed$conservation %in% c(1, 2)]
mature_mir_nonseed_conserved <- mature_mir_nonseed[mature_mir_nonseed$conservation %in% c(1, 2)]

## Only non-conserved
mature_mir_seed_nonconserved <- mature_mir_seed[mature_mir_seed$conservation %in% c(0)]
mature_mir_nonseed_nonconserved <- mature_mir_nonseed[mature_mir_nonseed$conservation %in% c(0)]

## Get 15bp flanks of primary miRNAs and remove any overlap with exons + 5bp to
## remove splice site influence
flk <- flank(primary_mir, width = 15)
flk_conserved <- flk[flk$ID %in% mature_mir_seed_conserved$Derives_from]
flk_nonconserved <- flk[flk$ID %in% mature_mir_seed_nonconserved$Derives_from]
flk_conserved <- setdiff(flk_conserved, exon_gr + 5, ignore.strand = TRUE)
flk_nonconserved <- setdiff(flk_nonconserved, exon_gr + 5, ignore.strand = TRUE)

## Export bed files
export_bed(mature_mir_seed, "mature_mir_seed.bed.gz", path = out_path)
export_bed(mature_mir_nonseed, "mature_mir_nonseed.bed.gz", path = out_path)
export_bed(mature_mir_seed_conserved, "mature_mir_seed_conserved.bed.gz", path = out_path)
export_bed(mature_mir_nonseed_conserved, "mature_mir_nonseed_conserved.bed.gz", path = out_path)
export_bed(mature_mir_seed_nonconserved, "mature_mir_seed_nonconserved.bed.gz", path = out_path)
export_bed(mature_mir_nonseed_nonconserved, "mature_mir_nonseed_nonconserved.bed.gz", path = out_path)
export_bed(flk_conserved, "flk_conserved.bed.gz", path = out_path)
export_bed(flk_nonconserved, "flk_nonconserved.bed.gz", path = out_path)
