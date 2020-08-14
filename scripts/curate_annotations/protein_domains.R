library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/protein_domains"
dir.create(out_path, FALSE, TRUE)

message("Loading protein domain annotations...")
## Load database and set filters
txdb <- load_txdb("GRCh38")
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter)

## Extract protein domain locations
prot <- proteins(txdb, columns = c("protein_id", "tx_id",
                           listColumns(txdb, "protein_domain")),
                 filter = filters)
prot <- prot[!is.na(prot$prot_dom_start) & !is.na(prot$prot_dom_end),]

## Subset to only pfam annotations
prot <- prot[prot$protein_domain_source == "pfam",]
prot <- prot[prot$interpro_accession != "",]

message("Filtering out rare domains...")
## Get the domains that cover at least 100,000 bp and are present at least 100 times
min_cover <- 1e5
min_count <- 100
domain_site_coverage <- integer(length(unique(prot$interpro_accession)))
names(domain_site_coverage) <- unique(prot$interpro_accession)
for(d in names(domain_site_coverage)) {
    domain_site_coverage[d] <-
        sum(with(prot[prot$interpro_accession == d,], prot_dom_end - prot_dom_start))
}

domain_count <- table(prot$interpro_accession)
common_domains <- intersect(names(domain_site_coverage >= min_cover),
                            names(domain_count[domain_count >= min_count]))

## Subset prot to only contain those common domains
common_prot <- prot[prot$interpro_accession %in% common_domains,]

txdb

## Convert protein domain locations to genome co-ordinates
message("Converting to genomic coordinates...")
prot_ir <- with(common_prot, IRanges(start = prot_dom_start, end = prot_dom_end))
mcols(prot_ir) <- with(common_prot,
                       DataFrame(protein_id, protein_domain_id, interpro_accession, protein_domain_source))
prot_gr <- GRangesList(proteinToGenome(prot_ir[1:1000], txdb, id = "protein_id"))
names(prot_gr) <- mcols(prot_ir[1:1000])$interpro_accession

## Filter out the bits with CDS issues
prot_all <- unlist(prot_gr, recursive = TRUE, use.names = TRUE)
prot_all <- prot_all[prot_all$cds_ok]

message("Writing bed files by domain...")
## Split by protein domain
prot_split <- split(prot_all, names(prot_all))
domains <- names(prot_split)
for(i in seq_along(domains)) {
    d <- domains[i]
    out_name <- paste0(d, ".bed.gz")    
    out_bed <- reduce(prot_split[[d]], ignore.strand = TRUE)
    export_bed(out_bed, file_name = out_name, path = out_path)
}

##############################
## library(data.table)
## protein_dat_dir <- "../../data/grch38/prot2hg"
## prot_feat_file <- file.path(protein_dat_dir, "protein_features_hg38.txt.gz")
## dir.create(protein_dat_dir, F, T)

## ## Download and save protein feature data
## ## pc_hg38_url <- "https://owncloud.cesnet.cz/index.php/s/uo97hRmf2DwBAqe/download"
## ## pc_dt <- fread(pc_hg38_url)
## ## fwrite(pc_dt, prot_feat_file)

## ## Load data
## prot_feat <- fread(prot_feat_file)

## ## Convert to GRanges
## with(prot_feat, GRanges(chrom, IRanges(chr_start, chr_end
