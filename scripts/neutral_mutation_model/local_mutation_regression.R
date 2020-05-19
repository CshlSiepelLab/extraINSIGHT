suppressWarnings(library(argparse))
#args = commandArgs(trailingOnly=TRUE)

parser <- ArgumentParser()
parser$add_argument("-b","--blocks", dest="blocks", type="character",
                    required=TRUE, help="bed file containing blocks for local normalization")
parser$add_argument("-m","--mutation-file", dest="mutation_file", type="character",
                    required=TRUE, help="file containing all mutations and their frequencies")
parser$add_argument("-n","--neutral", dest="neutral_file", type="character",
                    required=TRUE, help="neutral regions in bed format")
parser$add_argument("-g","--global-model", dest="global_model", type="character",
                    required=TRUE, help=".RData containing fitted GLM for global mutation model")
parser$add_argument("--flanking-size", dest="flanking_size", type="integer",
                    required=TRUE, help="length of flanking region to add to each block")
parser$add_argument("--min-coverage", dest="min_coverage", type="integer",
                    required=TRUE, help="minimum coverage filter for site inclusion")
parser$add_argument("--max-frequency", dest="max_frequency", type="double",
                    required=TRUE, help="maximum frequency filter for site inclusion")
parser$add_argument("--min-n-mutation", dest="min_n_mutation", type="integer",
                    required=TRUE, help="minimum number of neutral sites required to compute scaling factor for block")
parser$add_argument("-k","--context-frequency-table", dest="context_frequency_table", type="character",
                    required=TRUE, help="a table containing the frequencies of each k-mer context")
parser$add_argument("-c","--covariate_file", dest="covariate_file", type="character",
                    required=TRUE, help="a table containing the frequencies of each k-mer context")

args <- parser$parse_args()

suppressWarnings(library(data.table))
suppressWarnings(library(stringi))

# Function lifted from https://stackoverflow.com/questions/44763056/is-there-an-r-equivalent-of-pythons-string-format-function
strformat = function(str, vals) {
    vars = stringi::stri_match_all(str, regex = "\\{.*?\\}", vectorize_all = FALSE)[[1]][,1]
    x = str
    for (i in seq_along(names(vals))) {
        varName = names(vals)[i]
        varCode = paste0("{", varName, "}")
        x = stringi::stri_replace_all_fixed(x, varCode, vals[[varName]], vectorize_all = TRUE)
    }
    return(x)
}

compute_local_mutation_rates <- function(d, min_n_mutation, global_model) {
    status <- "@PASS"
    out <- data.table()
    if(nrow(d) == 0){
        status <- "@no_mutations_remaining"
    } else {
        mutation_count <- sum(d$mutation_flag== 1 & d$chrom_neutral != '.')
        ## Now check that there are enough sites for which there are mutations that are in
        ## neutral regions
        if (mutation_count >= min_n_mutation){        
            ## Columns to write out if scaling successful
            out_cols <- c(colnames(d)[1:6], "mutation_scaling_factor")
            ## Get the predicted mutation rate for each site based on the global model
            d[, predicted_mutation_rate := predict(global_model, newdata=d)]
            ## Subset out the neutral sites
            d_neut <- d[chrom_neutral != ".", ]
            ## Compute the scaling factor necessary to make the global model predictions fit the observed local neutral
            ## mutation rate as well as possible
            local_model <- glm(mutation_flag~predicted_mutation_rate, data=d_neut, family="binomial")
            d[,mutation_scaling_factor:= predict(local_model, newdata=d, type="response")]
            ## Write out the data
            out <- d[, ..out_cols]
        } else {
            status <- "@insufficient_neutral_mutations"
        }
        return(list(out = out, status = status))
    }
}

## load global model
load(args$global)
global_model <- model

## Load in blocks
blocks <- fread(args$blocks, header = FALSE)
blocks[, status := "@PENDING"]

for(i in seq_len(nrow(blocks))) {
    chr <- blocks[i, ]$V1
    bin_start <- blocks[i, ]$V2
    bin_end <- blocks[i, ]$V3
   
    # Add flanking regions to co-ordinates
    start <- max(0, bin_start - args$flanking_size) # if start is less than 0 make it 0
    end <- bin_end + args$flanking_size

    ## tabix queries assume a 1 based coordinate system (when not feeding in an input file with a *.bed or *.bed.gz suffix) 
    ## and our other co-ordinates are 0-based so we need to correct for that here so that the tabix queries are for the 
    ## correct co-ordinates. Outside the bed case tabix also assumes that the intervals are inclusive so there is no need
    ## to adjust the end co-ordinates.
    start <- start + 1
    
    ## Construct a command that extracts the mutations in the window +/- flanking regions and filters based on coverage and
    ## frequency. The remaining mutations are then merged with the genomic covariates ans subsetted again to contain only
    ## mutations in neutral regions. This subcommand returns a file with 13 columns:
    ## 1)chrom
    ## 2)start
    ## 3)end
    ## 4)context-alternate_allele
    ## 5)mutation_flag (0 = no mutation, 1 = mutation at that particular site)
    ## 6)average sequencing coverage
    ## 7)log(mutation probability) based on context from mutation_frequency_table.txt
    ## 8)GC_content
    ## 9)is_cpg_island
    ## 10)number of bases that this site overlaps with the the covariate bed file by
    ## 11)chrom (of overlapping neutral region, or "." if not in neutral region) 
    ## 12)start (of overlapping neutral region, or "-1" if not in neutral region)
    ## 13)end (of overlapping neutral region, or "-1" if not in neutral region)
    
    args_local <- c(args, list(chr = chr, start=start, end=end))
    str = paste("tabix {mutation_file} {chr}:{start}-{end} | perl filter_mutation.pl {min_coverage} {max_frequency} |",
        "perl add_context_mutation_rate.pl {context_frequency_table} | bedtools intersect -a - -b {covariate_file} -sorted -wo |",
        "cut -f1-7,11- | intersectBed -loj -a - -b {neutral_file} -sorted")
    cmd <- strformat(str, args_local)
    ## Run command and read into data.table 
    d <- fread(cmd = cmd, header = FALSE)
    ## Check that data read in as expected
    if(!exists("d")){
        stop("Failed to read in mutation data")
    }
    if(ncol(d) != 13){
        message("Malformed input (top row printed):")
        message(head(d, 1))
        stop(paste("Something has changed in data processing, input has",ncol(d),
                   "columns, should have 13. Has a covariate been added?"))
    }
    ## Rename columns
    setnames(d, colnames(d),
             c("chrom",
               "start",
               "end",
               "context",
               "mutation_flag",
               "coverage",
               "logit_mutation_rate",
               "gc_content",
               "cpg_island",
               "base_overlap",
               "chrom_neutral",
               "start_neutral",
               "end_neutral"
               )
             )
    ## Compute the adjusted mutations rate and get a status message per block
    results <- compute_local_mutation_rates(d, args$min_n_mutation, global_model)
    blocks[i, status := results$status]
    if (results$status == "@PASS") {
        fwrite(results$out, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }
}

fwrite(blocks, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
