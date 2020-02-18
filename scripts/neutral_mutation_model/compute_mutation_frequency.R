suppressPackageStartupMessages(library("argparse"))

## Setup parser
parser <- ArgumentParser()
parser$add_argument("-s", "--symmetrize-rates", action="store_true", default=FALSE, dest = "sym_rate",
                    help="Sum mutation counts from reverse-complementary mutation contexts with complementary mutations. [default %(default)s]")
parser$add_argument("input", nargs=1, help="input file of context incidence and mutation counts")

args <- parser$parse_args()

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))

## Read in data
mutation_counts = fread(paste("zcat",args$input))

if(args$sym_rate){
    ## Turn the mutation string into a key
    setkey(mutation_counts,"V1")

    ## Parse out context and alternate allele
    seq_matrix = with(mutation_counts,(str_split(V1,"-",simplify=T)))

    ## Get reverse complement of sequences
    reverse_context = as.character(reverseComplement(DNAStringSet(seq_matrix[,1])))
    reverse_mutation = as.character(reverseComplement(DNAStringSet(seq_matrix[,2])))

    ## Create a map of equivalent sequences 
    equivalence_map= data.table(mutation_counts$V1,paste0(reverse_context,"-",reverse_mutation))
    
    ## Add counts of sequences that are equivalent
    symm_mut_rate = cbind(equivalence_map$V1,mutation_counts[equivalence_map$V1,.(V2,V3)] +
        mutation_counts[equivalence_map$V2,.(V2,V3)])
    ## Compute symmetrized mutation rate
    symm_mut_rate[,V4:=V3/V2]
    ## Create output table
    mutation_rate_out = symm_mut_rate[,.(V1,V4)]
} else {
    ## Compute mutation rate
    mutation_rate_out = mutation_counts[,.(V1,V3/V2)]
}

## Make sure output is sorted by mutation string
mutation_rate_out = mutation_rate_out[order(V1),]

## Write to stdout
fwrite(mutation_rate_out, file = "", col.names = F,sep = "\t")
