library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
source("annotation_lib.R")

## Load in annotation database
args = list()
args$txdb = "grch38.sqlite"
args$out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/reactome"
txdb <- loadDb(args$txdb)
txdb <- AnnotationDbi::loadDb(args$txdb)
txdb <- keepSeqlevels(txdb, as.character(1:22))
dir.create(args$out_dir, FALSE, TRUE)

## Get coding genes
tx <- transcripts(txdb, columns=c("GENEID", "TXTYPE", "TXNAME"))
coding_tx <- unlist(tx[tx$TXTYPE == "protein_coding"]$TXNAME)
coding_gene <- unique(unlist(tx[tx$TXTYPE == "protein_coding"]$GENEID))
tx_gene_map <- with(tx, data.table(GENEID = unlist(GENEID), TXNAME, TXTYPE, key = "GENEID"))

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
       file = file.path(args$out, "id_pathway_mapping.txt.gz"))
gene_to_pathway_final <- gene_to_pathway[.(final_terms)]

## Export bed files for each top level reactome annotation
## CDS
for (i in final_terms) {
    gene_subset <- gene_to_pathway[reactome_id == i]$source_db_id
    gr_out <- cds(txdb, columns = "GENEID", filter = list(GENEID = gene_subset))
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(gr_out, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_cds.bed.gz"), args$out_dir)
}

## 5' UTR for protein coding transcripts only
five_gr <- fiveUTRsByTranscript(txdb, use.names = TRUE)
five_gr <- unlist(five_gr[intersect(coding_tx, names(five_gr))])
cds_gr <- cds(txdb) # use this to exclude regions that are ever coding
for (i in final_terms) {
    gene_subset <- gene_to_pathway[reactome_id == i]$source_db_id
    tx_subset <- tx_gene_map[gene_subset][TXTYPE == "protein_coding"]$TXNAME
    gr_out <- five_gr[intersect(tx_subset,names(five_gr))]
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(gr_out, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_fiveutr.bed.gz"), args$out_dir)
}

## 3' UTR for protein coding transcripts only
three_gr <- threeUTRsByTranscript(txdb, use.names = TRUE)
three_gr <- unlist(three_gr[intersect(coding_tx, names(three_gr))])
cds_gr <- cds(txdb) # use this to exclude regions that are ever coding
for (i in final_terms) {
    gene_subset <- gene_to_pathway[reactome_id == i]$source_db_id
    tx_subset <- tx_gene_map[gene_subset][TXTYPE == "protein_coding"]$TXNAME
    gr_out <- three_gr[intersect(tx_subset,names(three_gr))]
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(gr_out, ignore.strand = TRUE)
    export_bed(gr_out, paste0(i, "_threeutr.bed.gz"), args$out_dir)
}
