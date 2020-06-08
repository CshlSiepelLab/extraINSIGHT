## Calculate annotations based on gene expression per tissue
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

## Ensembl release 99 = Gencode v33
## txdb <- makeTxDbFromEnsembl(organism="Homo sapiens", release = 99)
## saveDb(txdb, file="grch38.sqlite")

args = list()
args$txdb = "grch38.sqlite"
args$out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/top_5_expressed"

txdb <- loadDb(args$txdb)

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

txdb <- AnnotationDbi::loadDb(args$txdb)
cds_gr <- sort(cds(txdb,columns=c("GENEID")))
cds_gr$GENEID <- unlist(cds_gr$GENEID)

## Get paths and controlled vocalbulary for gtex data 
gtex_data_path <- "../../data/grch38/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

# Read in gtex data
tpm <- fread(gtex_data_path, header = TRUE)
tpm_melted <- melt(tpm, id.vars = c("Name", "Description"),
                   variable.name = "Tissue")
# Compute percentile of each gene among non-zero genes
tpm_melted[value > 0,  quantile := order(value) / length(value), by = "Tissue"]
setkey(tpm_melted, "Tissue")
tpm_melted[,Name := gsub("\\.\\d+", "", Name)]
dir.create(args$out_dir, FALSE, TRUE)

for(tis in as.character(unique(tpm_melted$Tissue))) {
    label <- gsub("-","",tis)
    label <- gsub("\\s+","_",label)
    label <- gsub("\\(|\\)","",label)
    tmp <- cds_gr[cds_gr$GENEID %in% tpm_melted[tis][quantile >= 0.95]$Name]
    tmp <- keepSeqlevels(tmp, 1:22, pruning.mode = "coarse")
    export_bed(tmp, paste0(label,"_cds.bed.gz"), path = args$out_dir)
}
