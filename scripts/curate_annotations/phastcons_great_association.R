library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

basal_domains <- function(gr, upstream, downstream) {
    basal <- promoters(gr, upstream, downstream)
    strand(basal) <- "*"
    return(basal)
}

regulatory_domains <- function(gr, upstream = 5e3, downstream = 1e3, distal_max = 1e6) {
    ## Get basal domains
    bd <- basal_domains(gr, upstream, downstream)
    bd <- sort(sortSeqlevels(bd))

    ## Take one bp off to each side of basal region
    bp_left <- GRanges(seqnames(bd), IRanges(start(bd) - 1, start(bd) - 1))
    bp_right <- GRanges(seqnames(bd), IRanges(end(bd) + 1, end(bd) + 1))

    ## Check if bp intersects with another basal domain
    do_not_extend_left <-  queryHits(findOverlaps(bp_left, bd))
    do_not_extend_right <- queryHits(findOverlaps(bp_right, bd))

    ## Compute the start/end of next/previous basal domain for each gene with fixes for first and
    ## last genes on each chromosomes
    previous_basal_end <- end(bd)[follow(bp_left, bd)]
    previous_basal_end[is.na(previous_basal_end)] <- 0
    next_basal_start <- start(bd)[precede(bp_right, bd)]
    next_basal_start[is.na(next_basal_start)] <-
        seqlengths(bd)[decode(seqnames(bp_right[is.na(next_basal_start)]))]
    
    ## Compute extensions with upper limit on length
    left_ext <- pmin(start(bd) - previous_basal_end, distal_max) - 1
    right_ext <- pmin(next_basal_start - end(bd), distal_max) - 1
    left_ext[do_not_extend_left] <- 0
    right_ext[do_not_extend_right] <- 0

    ## Create regulatory domains
    r_dom <- GRanges(seqnames(bd), IRanges(start(bd) - left_ext,
                                           end(bd) + right_ext))
    mcols(r_dom) <- mcols(bd)
    return(r_dom)
}

out_dir <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/phastcons_functional_annotations"
unlink(out_dir, T)
dir.create(out_dir, FALSE, TRUE)

## Load in genes
txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)
gene <- genes(txdb, filter = filters)
coding_gene <- gene$gene_id

##
## Annotate genes with reactome top level pathways
##

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


##
## Associate genes and phastCons elements
##

## Load in phastCons elements
phastcons <- "../../data/grch38/annotation_bed/phastcons_elements/phastcons_elements_100way.bed.gz"
bed <- with(fread(phastcons), GRanges(V1, IRanges(V2, V3-1)))
seqlevelsStyle(bed) <- "ensembl"

## Create regulatory domains and overlap regions with them
reg_domains <- regulatory_domains(gene)


## Subset out each reactome catagory and the accompanying genes and find the
## phastcons elements that overlap with their regulatory domains
setkey(gene_to_pathway_final, "reactome_id")
for(i in unique(gene_to_pathway_final$reactome_id)) {
    gene_subset <- gene_to_pathway_final[i]$source_db_id
    associations <- findOverlaps(bed, reg_domains[reg_domains$gene_id %in% gene_subset])
    phastcons_sub <- unique(queryHits(associations))
    i <- gsub("\\s+", "_", i)
    gr_out <- reduce(bed[phastcons_sub])
    export_bed(gr_out, paste0(i, "_phastcons.bed.gz"), out_dir)
}
