library(GenomicRanges)
library(readxl)
library(GenomeInfoDb)
source("annotation_lib.R")

out_path <- "~/Projects/extraINSIGHT/data/grch37/annotation_bed/har_doan"
unlink(out_path, T)
dir.create(out_path, FALSE, TRUE)

## Download supplementary table 1 from doan et al. 2016
## <https://www.cell.com/fulltext/S0092-8674(16)31169-2#supplementaryMaterial>
## Table S1. Combined HARs from Five Previous Studies in hg19 Genome Build
## Annotated with ChromHMM Predictions, refSeq Genes, and Target Genes from
## Existing ChIA-PET and HiC Data, Related to Figure 1
har_url <- "https://www.cell.com/cms/10.1016/j.cell.2016.08.071/attachment/9e010545-053a-4166-9b0d-7cb8c8f90ca7/mmc2.xlsx"
tmp <- tempfile()
download.file(har_url, tmp)

## Read in excel sheet
har <- read_excel(path = tmp, sheet = 1)

## Convert to GRanges
har_gr <- with(har, GRanges(Chr, IRanges(Start, End - 1)))
har_gr <- keepSeqlevels(har_gr, paste0("chr",as.character(1:22)), pruning.mode = "coarse")

## Export
export_bed(reduce(har_gr), "har_doan.bed.gz", path = out_path)
