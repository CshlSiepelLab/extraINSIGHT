#! /usr/bin/env Rscript

library(Biostrings)
library(R.utils)
library(data.table)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ensembldb)
source("annotation_lib.R")
source("load_txdb.R")

## Function that pre-calculates table that give the degree of degeneracy at
## each site in the codon (ranges from 1(no substitutions) - 4(any substitution))
calculate_degeracy_degree <- function(genetic_code = GENETIC_CODE){
    genetic_code <- as.list(genetic_code)
    degeneracy_degree <- lapply(genetic_code, function(x) integer(3))
    alp <- c("A", "C", "G", "T")
    for(g in names(genetic_code)){
        aa <- genetic_code[[g]]
        for(p in seq_len(3)){
            for(a in alp){
                tmp <- g
                substr(tmp, p, p) <- a
                degeneracy_degree[[g]][p] <- degeneracy_degree[[g]][p] +
                    as.integer(genetic_code[[tmp]] == aa)
            }
        }
    }
    return(degeneracy_degree)
}


## Function for coverting an RLE to a granges
## Adapted from:
## https://stackoverflow.com/questions/39222913/efficiently-construct-granges-iranges-from-rle-vector
rle_to_GRanges <- function(seqname, r, strand = NULL, name = NULL){
    gr <- GRanges(rep(seqname, length.out = length(runValue(r))),
                  IRanges(cumsum(c(0,runLength(r)[-nrun(r)])),
                          width=runLength(r), names = name),
                  strand = strand,
                  score = runValue(r))
    return(ir)    
}

rle_to_IRanges <- function(r){
    ir <-  IRanges(cumsum(c(0,runLength(r)[-nrun(r)])),
                   width=runLength(r),
                   score = runValue(r))
    return(ir)
}


## Vectorized seq function
seqv <- Vectorize(seq.default, vectorize.args = c("from", "to"))

## Create a function for mapping from CDS coordinates to GenomeCoordinates
## ir = IRanges for a single CDS (start at 1)
## transcript_id = ENST id
## edb = ensembldb
cds_to_genome_coordinate <- function(ir, transcript_id, edb) {
    tx <- cdsBy(edb, filter = TxIdFilter(transcript_id))[[1]]
    sqnm <- runValue(seqnames(tx))
    strnd <- runValue(strand(tx))
    ## Get per bp indicies of cds
    per_bp_cds <- Reduce("c",seqv(start(ir), end(ir)))
    ## Create a vector of values where the index of the entry maps to the
    ## genomic coordinate
    if(strnd == "-") {
        arr <- Reduce("c", seqv(from = end(tx), to = start(tx), by = -1))[per_bp_cds]
    } else {
        arr <- Reduce("c", seqv(from = start(tx), to = end(tx), by = 1))[per_bp_cds]
    }
    ## Create GRanges of converted genomic coordinates
    gr <- sort(GRanges(sqnm, IRanges(arr, width = 1), strand = strnd,
                       score = rep(mcols(ir)$score, times = width(ir))))
    ## Split gr by score so that adjacent within score ranges can be merged
    gr <- split(gr, score(gr))
    gr <- unlist(reduce(gr, min.gapwidth=1L), use.names = TRUE)
    score(gr) <- names(gr)
    gr <- sort(gr)
    names(gr) <- rep(transcript_id, length(gr))
    return(gr)
}

## Calculate degeneracy table
deg_table <- calculate_degeracy_degree()

## Out path
out_path <- "~/Projects/extraINSIGHT/data/grch38/annotation_bed/4d_sites"
dir.create(out_path, FALSE, TRUE)

## Load Ensdb
txdb <- load_txdb("GRCh38")

## Filters

## TSL Categories
## The following categories are assigned to each of the evaluated annotations:
## tsl1  all splice junctions of the transcript are supported by at least one non-suspect mRNA
## tsl2  the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
## tsl3  the only support is from a single EST
## tsl4  the best supporting EST is flagged as suspect
## tsl5  no single transcript supports the model structure
txs_filter <- TxSupportLevelFilter(1)
seqnames_filter <- SeqNameFilter(as.character(1:22))
txbiotype_filter <- TxBiotypeFilter("protein_coding")
filters <- AnnotationFilterList(seqnames_filter, txbiotype_filter, txs_filter)

## Get list of genes to subsample
set.seed(3409)
subsample <- 2e3
genes <- genes(txdb, filter = filters)
sub_idx <- sample(seq_along(genes), subsample)
gene_id_sub <- genes[sub_idx]$gene_id
gene_filter <- GeneIdFilter(gene_id_sub)
filters <- c(AnnotationFilterList(gene_filter), filters)

## Extract a variety of regions
## NOTE: this retrieves exons on the order of transcription (exon rank)
## so in order to get them in the correct order for retrieving sequences that can be cleanly
## pasted together (and reverse complemented when necessary), ranges are resorted.
## Also extract the protein sequence so that it can serve as a QC step to remove sequences
## that did not translate correctly
cds_gr <- unlist(sort(cdsBy(txdb, filter = filters, by = "tx",
                            columns = c("protein_sequence", "protein_id"))))
## De-duplicate protein sequences such that you get one per transcript
ebl_protein_seq <- cds_gr$protein_sequence[which(!duplicated(names(cds_gr)))]
names(ebl_protein_seq) <- names(cds_gr)[which(!duplicated(names(cds_gr)))]

## Get twoBit file
fa <- getGenomeTwoBitFile(txdb)

## Extract CDS per transcript
cds_split_seq <- import(fa, which = cds_gr)

## Concatonate CDS to a single sequence
cds_joint <- DNAStringSet(
    lapply(split(cds_split_seq, names(cds_split_seq)), function(x)
           DNAString(paste(as.character(x), collapse = "")))
    )

## Extract strand and reverse complement where "-"
cds_strand <- unlist(runValue(RleList(split(strand(cds_gr), names(cds_split_seq)))))
cds_seqname <- unlist(runValue(RleList(split(seqnames(cds_gr), names(cds_split_seq)))))
cds_joint[cds_strand == "-"] <- reverseComplement(cds_joint[cds_strand == "-"])
aa_joint <- translate(cds_joint)

## Check that infered protein sequence matches the one in the ensembl database
## and remove those that don't (often b/c of a 5' truncation or something)
keep_tx <- names(aa_joint)[gsub("\\*", "",as.character(aa_joint)) == ebl_protein_seq]
cds_joint <- cds_joint[keep_tx]

## Split sequences into codons and replace each site with its degree of degeneracy
seq_deg <- IRangesList(lapply(cds_joint, function(x){
    r <- Rle(unlist(deg_table[stri_sub(x, seq(1, stri_length(x),by=3), length=3)], use.names = F))
    return(rle_to_IRanges(r))
}))

## Map to genome coordinates
genome_coord_deg <- GRangesList(BiocGenerics::mapply(function(x, name){
    x <- shift(x , 1)
    return(cds_to_genome_coordinate(x,name, txdb))
}, x = seq_deg, name = as.list(names(seq_deg))))

all_deg_gr <- unlist(genome_coord_deg, use.names = FALSE)
score(all_deg_gr) <- as.integer(score(all_deg_gr))

## Get 4-fold degenerate sites
fourfold <- all_deg_gr[score(all_deg_gr) == 4]
## And non-degenrate sites
oned <- all_deg_gr[score(all_deg_gr) == 1]

## Make sure that there are no overlaps
fourfold <- setdiff(fourfold, oned)
oned <- setdiff(oned, fourfold)

export_bed(fourfold, "fourfold_degenerate_sites.bed.gz", out_path)
export_bed(oned, "non_degenerate_sites.bed.gz", out_path)

## Some ideas for how to batch qery the ensem

## library(httr)
## library(jsonlite)
## library(xml2)

## fetch_endpoint_POST <- function(server, endpoint, data, content_type='application/json'){
##     r <- POST(paste0(server,"/" ,endpoint),
##               content_type("application/json"),
##               accept("application/json"),
##               body = data)

##     stop_for_status(r)

##     if (content_type == 'application/json'){
##         return (fromJSON(content(r, "text", encoding = "UTF-8")))
##     } else {
##         return (content(r, "text", encoding = "UTF-8"))
##     }
## }


## fetch_endpoint_POST(server, endpoint, body_values)

## server <- "https://rest.ensembl.org"
## endpoint <- "map/cds/"

## start <- start(x)
## end <- end(x)
## gene <- rep("ENST00000001146", length(start))
## region <- paste0(start, "..", end)
## body_values <- toJSON(list(ids=gene, regions = region))





 
## # use this if you get a simple nested list back, otherwise inspect its structure
## # head(data.frame(t(sapply(content(r),c))))
## head(fromJSON(toJSON(content(r))))





