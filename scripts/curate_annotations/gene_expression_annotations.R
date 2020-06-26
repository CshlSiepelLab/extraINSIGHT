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

args <- list()
args$out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/top_5_expressed"

## Setup txdb and filters
txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
tx_biotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, tx_biotype_filter)

## Load in coding genes and 
coding_genes <- genes(txdb, filter = filters)
cds_gr <- unlist(cdsBy(txdb, filter = filters, by = "gene"), use.names = TRUE)

##################################################
## Annotate highly expressed genes in each tissue
##################################################

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
dir.create(args$out_dir, FALSE, TRUE)

for(tis in as.character(unique(tpm_melted$Tissue))) {
    label <- gsub("-","",tis)
    label <- gsub("\\s+","_",label)
    label <- gsub("\\(|\\)","",label)
    tmp <- reduce(cds_gr[intersect(names(cds_gr), tpm_melted[tis][quantile >= 0.95]$Name)])
    export_bed(tmp, paste0(label,"_cds.bed.gz"), path = args$out_dir)
}

###################
## Pleiotropy
###################
args$out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/expression_pleiotropy"
plots_dir <- file.path(args$out_dir, "plots")
dir.create(plots_dir, FALSE, TRUE)

expression_cutoff <- 0.1
pleio <- tpm_melted[, .(num_tissues_expressed = sum(value >= expression_cutoff)), by = "Name"]
q <- quantile(pleio$num_tissues_expressed, seq(0,1,0.05))
pleio[num_tissues_expressed == 0, category := "unexpressed"]
pleio[num_tissues_expressed <= 5 & num_tissues_expressed > 0, category := "tissue specific"]
pleio[num_tissues_expressed >= 50, category := "broadly expressed"]
pleio[is.na(category), category := "intermediate"]

# Plot out distribution of number of tissues expressed
g <- ggplot(pleio, aes(x = num_tissues_expressed, fill = category))+
    geom_histogram(bins = 55)+
    theme_cowplot()+
    theme(legend.position=c(0.6,0.8), legend.box.background = element_rect(colour = "black"))+
    scale_fill_viridis_d()+
    scale_x_continuous(breaks = pretty(0:54, n = 10))+
    xlab(paste0("Tissues expressed (TPM >= ", expression_cutoff, ")"))+
    ylab("Frequency")+
    guides(fill = guide_legend(title = "Category"))
ggsave("number_tissues_expressed.pdf", path = plots_dir, plot = g, height = 5, width = 6)

## CDS
for (i in unique(pleio$category)) {
    gene_subset <- pleio[category == i]$Name
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    gr_out <- reduce(sort(unlist(cdsBy(txdb, filter = iter_filters))))
    i <- gsub("\\s+", "_", i)
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), args$out_dir)
}

## 5' UTR for protein coding transcripts only
cds_gr <- unlist(cdsBy(txdb, filter = filters)) # use this to exclude regions that are ever coding
for (i in unique(pleio$category)) {
    gene_subset <- pleio[category == i]$Name
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    five_gr <- unlist(fiveUTRsByTranscript(txdb, filter = iter_filters))
    five_gr <- reduce(setdiff(five_gr, cds_gr))
    i <- gsub("\\s+", "_", i)
    export_bed(gr_out, paste0(i, "_fiveutr.bed.gz"), args$out_dir)
}
                 
## 3' UTR for protein coding transcripts only
cds_gr <- unlist(cdsBy(txdb, filter = filters)) # use this to exclude regions that are ever coding
for (i in unique(pleio$category)) {
    gene_subset <- pleio[category == i]$Name
    gene_id_filter <- GeneIdFilter(gene_subset)
    iter_filters <- c(filters, AnnotationFilterList(gene_id_filter))
    three_gr <- unlist(threeUTRsByTranscript(txdb, filter = iter_filters))
    three_gr <- reduce(setdiff(three_gr, cds_gr))
    i <- gsub("\\s+", "_", i)
    export_bed(gr_out, paste0(i, "_threeutr.bed.gz"), args$out_dir)
}
