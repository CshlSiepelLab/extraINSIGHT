## Calculate annotations based on gene expression per tissue
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ggplot2)
library(cowplot)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

## Setup txdb and filters
txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
tx_biotype_filter <- TxBiotypeFilter("protein_coding")
gene_biotype_filter <- GeneBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, tx_biotype_filter)

## Load in coding genes and 
coding_genes <- genes(txdb, filter = filters)
cds_gr <- unlist(cdsBy(txdb, filter = filters, by = "gene"), use.names = TRUE)

##################################################
## Annotate highly expressed genes in each tissue
##################################################
out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/top_5_expressed"
unlink(out_dir, TRUE)
## Get median gene expression
gtex_data_path <- "../../data/grch38/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

## Read in gtex data and filter out all non-coding genes
tpm <- fread(gtex_data_path, header = TRUE)
tpm_melted <- melt(tpm, id.vars = c("Name", "Description"),
                   variable.name = "Tissue")
tpm_melted[,Name := gsub("\\.\\d+", "", Name)]
setkey(tpm_melted, "Name")
tpm_melted <- tpm_melted[intersect(coding_genes$gene_id, tpm_melted$Name)]
## Compute percentile of each gene among non-zero genes
tpm_melted[value > 0,  quantile := order(value) / length(value), by = "Tissue"]
setkey(tpm_melted, "Tissue")
dir.create(out_dir, FALSE, TRUE)

for(tis in as.character(unique(tpm_melted$Tissue))) {
    label <- gsub("-","",tis)
    label <- gsub("\\s+","_",label)
    label <- gsub("\\(|\\)","",label)
    tmp <- reduce(cds_gr[intersect(names(cds_gr), tpm_melted[tis][quantile >= 0.95]$Name)])
    export_bed(tmp, paste0(label,"_cds.bed.gz"), path = out_dir)
}

################################################################
## Pleiotropy - number of tissues expressed at above median TPM
################################################################
out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/expression_pleiotropy"
unlink(out_dir, TRUE) 
plots_dir <- file.path(out_dir, "plots")
dir.create(plots_dir, FALSE, TRUE)
dir.create(out_dir, FALSE, TRUE)
## Get median gene expression
gtex_data_path <- "../../data/grch38/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

## Read in gtex data and filter out all non-coding genes
tpm <- fread(gtex_data_path, header = TRUE)
tpm_melted <- melt(tpm, id.vars = c("Name", "Description"),
                   variable.name = "Tissue")
tpm_melted[,Name := gsub("\\.\\d+", "", Name)]
setkey(tpm_melted, "Name")
tpm_melted <- tpm_melted[intersect(coding_genes$gene_id, tpm_melted$Name)]

over_median <- tpm_melted[, .(Name = Name, over=value >= median(value)), by = "Tissue"]
tissues_expressed <- over_median[, .(num_tissues_expressed = sum(over)), by = "Name"]

tissues_expressed[num_tissues_expressed == 0, category := "unexpressed"]
tissues_expressed[num_tissues_expressed <= 5 & num_tissues_expressed > 0, category := "tissue specific"]
tissues_expressed[num_tissues_expressed >= 50, category := "broadly expressed"]
tissues_expressed[is.na(category), category := "intermediate"]

## CDS
for (i in unique(tissues_expressed$category)) {
    gene_subset <- tissues_expressed[category == i]$Name
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    gr_out <- reduce(sort(unlist(cdsBy(txdb, filter = iter_filters))))
    i <- gsub("\\s+", "_", i)
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), out_dir)
}

## # Plot out distribution of number of tissues expressed
g <- ggplot(tissues_expressed, aes(x = num_tissues_expressed, fill = category))+
    geom_histogram(bins = 55)+
    theme_cowplot()+
    theme(legend.position=c(0.6,0.8), legend.box.background = element_rect(colour = "black"))+
    scale_fill_viridis_d()+
    scale_x_continuous(breaks = pretty(0:54, n = 10))+
    xlab(paste0("Tissues expressed (TPM >= tissue_median)"))+
    ylab("Frequency")+
    guides(fill = guide_legend(title = "Category"))
ggsave("number_tissues_expressed.pdf", path = plots_dir, plot = g, height = 5, width = 6)

###################
## Tissue specific expression
###################
out_dir <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/tissue_specific_expression"
unlink(out_dir, T)
txdb_37 <- load_txdb("GRCh37")
dir.create(out_dir, F, T)
## Load tissue specificity scores from supplemental table 1 of yang (2018) -bioarxiv
## https://www.biorxiv.org/content/10.1101/311563v1.full#F1
## Elements are only in file when they have TS >= 3 which means that they are
## "tissue specific" for that tissue. That is roughly equivalent to having 2^3 times
## more expression in that tissue than over a weighted average of all the other tissues.
## Genes can be "tissue specific" for multiple tissues. They used GTEx 6p so we do too for
## this labeling.
ts_scores <- fread("../../data/grch38/gtex/ts_score_stab2_yang2018.csv.gz")
ts_scores[, Name:= gsub("\\.\\d+$", "", Name)]
ts_scores <- ts_scores[type == "protein_coding"]

## Get list of genes used in analysis for v6 of GTEx
gtex_v6p_url <- paste0("https://storage.googleapis.com/gtex_analysis_v6/reference/",
                       "gencode.v19.genes.patched_contigs.gtf.gz")
v6p_annotations <- as.data.table(rtracklayer::readGFF(gtex_v6p_url))

## Print tissue specific gene CDS for each tissue
for (i in unique(ts_scores$tissue)) {
    gene_subset <- ts_scores[tissue == i]$Name
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    gr_out <- reduce(sort(unlist(cdsBy(txdb_37, filter = iter_filters))))
    i <- tolower(gsub("\\s+|-", "_", i))
    i <- gsub("_+","_", i)
    i <- gsub("\\(|\\)", "", i)
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), out_dir)
}

###################
##  Pleiotropy based on Tissue specific expression
###################
out_dir <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/tissue_group_exclusivity"
unlink(out_dir, T)
dir.create(out_dir, F, T)
## Supplementary Table 8. A table of TS_Score of all genes in 49 tissues.
## TS_Score is the tissue specificity score. Group.cnt is the  number of
## tissue groups the gene shows enriched expression (i.e. TS_Score>3).
## TPM: gene expression value in transcripts per million.
##
## Now split protein coding genes into tissue exclusive, tissue specific,
## tissue general (low expression), tissue general (high expression) based
## on the number of tissue groups they are expressed in
pc_genes <- genes(txdb_37, filter = c(filters, AnnotationFilterList(gene_biotype_filter)))
all_ts <- fread("../../data/grch38/gtex/all_genes_ts_score_yang2018.txt.gz")
all_ts[, Name:= gsub("\\.\\d+$", "", Name)]
all_ts <- all_ts[Name %in% names(pc_genes)]

## Set types
tissue_group_specific <- unique(all_ts[groups.cnt > 1]$Name)
tissue_group_exclusive <- unique(all_ts[groups.cnt == 1]$Name)
median_tpm <- all_ts[!Name %in% c(tissue_group_specific, tissue_group_exclusive), median(tpm),by = "Name"]
## split on median of expressed genes
cutoff <- median(median_tpm$V1[median_tpm$V1> 0])
general_lo <- median_tpm[V1 > 0 & V1 <= cutoff]$Name
general_hi <- median_tpm[V1 > cutoff]$Name

cata <- list(tissue_group_specific, tissue_group_exclusive, general_hi, general_lo)
names(cata) <- c("tissue_group_specific", "tissue_group_exclusive", "general_hi", "general_lo")

## Print tissue specific gene CDS for each tissue
for (i in names(cata)) {
    gene_subset <- cata[[i]]
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    gr_out <- reduce(sort(unlist(cdsBy(txdb_37, filter = iter_filters))))
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), out_dir)
}

