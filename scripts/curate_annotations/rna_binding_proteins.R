library(data.table)
library(stringr)
library(GenomicRanges)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

## From postar2
## http://lulab.life.tsinghua.edu.cn/postar2
rnabp <- fread("../../data/grch38/postar2_rnabp/human_RBP_binding_sites.txt.gz")
cols <-c("chrom",
         "peak_start",
         "peak_end",
         "name",
         "0",
         "strand",
         "rbp_name",
         "clip-tech_peak-caller",
         "cell_line_tissue",
         "data_accession",
         "score")
setnames(rnabp, colnames(rnabp), cols)
rnabp[, peak_caller := str_match(`clip-tech_peak-caller`, ",(.*)")[,2]]

## Make into GRanges, minus one for conversion [,) --> [,]
rnabp_gr <- with(rnabp, GRanges(chrom, IRanges(peak_start, peak_end - 1, names = rbp_name),
                                strand = strand, cell_line_tissue, peak_caller, score))

############################################
## Export RNA binding protein binding sites
############################################
out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/rna_binding_proteins"
dir.create(out_path, FALSE, TRUE)


rnabp_gr_split <- split(rnabp_gr, names(rnabp_gr))
## Export as per RNA BP bed files
rnabp_names <- unique(names(rnabp_gr_split))
for(i in seq_along(rnabp_names)) {
    out_name <- paste0(rnabp_names[i], ".bed.gz")
    export_bed(rnabp_gr_split[[rnabp_names[i]]], out_name, out_path)
}
