library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
source("annotation_lib.R")
source("load_txdb.R")

## Load in annotation database
out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/reactome"
dir.create(out_dir, FALSE, TRUE)

txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
genebiotype_filter <- GeneBiotypeFilter("protein_coding")
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, genebiotype_filter)

## Get coding genes
genes <- genes(txdb, filter = filters)
coding_gene <- genes$gene_id

## Get pathway heirarchy
pathway_hierarchy_path <- "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
pathway_heirarchy <- fread(pathway_hierarchy_path, header = FALSE)
setnames(pathway_heirarchy, colnames(pathway_heirarchy), c("parent", "child"))
## Extract top level identifiers
top_level_terms <- setdiff(unique(pathway_heirarchy$parent), unique(pathway_heirarchy$child))

## Download terms to gene mapping
gene_to_pathway_path <- "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
gene_to_pathway <- fread(gene_to_pathway_path, header = FALSE,
                         col.names = c("source_db_id", "reactome_id", "url",
                             "pathway", "evidence", "species")
                         )[species == "Homo sapiens"]
gene_to_pathway <- gene_to_pathway[,.(source_db_id, reactome_id, pathway)]

## Keep only coding genes (only gene level)
setkeyv(gene_to_pathway, c("source_db_id"))
gene_to_pathway <- gene_to_pathway[.(coding_gene)]
gene_to_pathway <- gene_to_pathway[!is.na(reactome_id)]

## Remove duplicate entries
setkeyv(gene_to_pathway, c("reactome_id", "source_db_id"))
gene_to_pathway <- unique(gene_to_pathway)

## Subset to human only top level term entries
human_top_level_terms <- intersect(unique(gene_to_pathway$reactome_id), top_level_terms)
gene_to_pathway <- gene_to_pathway[.(human_top_level_terms)]

## Count genes per pathways and remove ones with less than 100 genes
gene_per_pathway <- gene_to_pathway[, .(
    reactome_id = reactome_id[1],
    gene_members = length(unique(source_db_id))), by = "pathway"]
final_terms <- gene_per_pathway[gene_members >= 100]$reactome_id
fwrite(gene_per_pathway[gene_members >= 100],
       file = file.path(out_dir, "id_pathway_mapping.txt.gz"))
gene_to_pathway_final <- gene_to_pathway[.(final_terms)]

## Export bed files for each top level reactome annotation
## CDS
for (i in final_terms) {
    gene_filter <- AnnotationFilterList(GeneIdFilter(gene_to_pathway[reactome_id == i]$source_db_id))
    gr_out <- reduce(unlist(cdsBy(txdb, filter = c(gene_filter, filters))))
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(gr_out, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), out_dir)
}

## 5' UTR for protein coding transcripts only
cds_gr <- unlist(reduce(cdsBy(txdb, filter = filters))) # use this to exclude regions that are ever coding
for (i in final_terms) {
    gene_filter <- AnnotationFilterList(GeneIdFilter(gene_to_pathway[reactome_id == i]$source_db_id))
    final_filter <- c(gene_filter, filters, AnnotationFilterList(txbiotype_filter))
    five_gr <- unlist(fiveUTRsByTranscript(txdb, filter = final_filter))
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(five_gr, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_fiveutr.bed.gz"), out_dir)
}

## 3' UTR for protein coding transcripts only
cds_gr <- unlist(reduce(cdsBy(txdb, filter = filters))) # use this to exclude regions that are ever coding
for (i in final_terms) {
    gene_filter <- AnnotationFilterList(GeneIdFilter(gene_to_pathway[reactome_id == i]$source_db_id))
    final_filter <- c(gene_filter, filters, AnnotationFilterList(txbiotype_filter))
    three_gr <- unlist(threeUTRsByTranscript(txdb, filter = final_filter))
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(three_gr, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_threeutr.bed.gz"), out_dir)
}
