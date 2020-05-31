suppressWarnings(library(argparse))

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
parser$add_argument("-o","--output", dest="output", type="character",
                    required=TRUE, help="path to write finalized mutation rates to")
parser$add_argument("-l","--log", dest="log", type="character",
                    required=TRUE, help="path to write block status to")
args <- parser$parse_args()

## tst <- "-b /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/jobs/window_chunk.bed.00042 -m /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/gnomad_all_potential_mutation.bed.gz -n /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/data/grch38/neutral_regions/grch38_neutral_region.bed -g /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/glm_slim.RData --flanking-size 25000 --min-coverage 20 --max-frequency 0.001 --min-n-mutation 200 -k /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/mutation_frequency_table.txt -c /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/data/grch38/covariates/unified_covariate_file.bed.gz -l /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/jobs_output/block_status.00250.bed -o /mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/jobs_output/mutation_rates.00250.bed.gz"
## tst <- unlist(strsplit(tst, " "))
## args <- parser$parse_args(tst)

suppressWarnings(library(data.table))
suppressWarnings(library(stringi))

# Function altered from https://stackoverflow.com/questions/7014081/capture-both-exit-status-and-output-from-a-system-call-in-r/17090601
robust_system <- function (cmd, retry_wait = 0, restart_times = 0) {
    retry <- TRUE
    while (retry) {
        tmp <- tryCatch({
            retval <- list()
            stderrFile <- tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
            stdoutFile <- tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
            retval$exitStatus <- system(paste0(cmd, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
            if (length(retval$stderr) > 0) {
                message("system call error:", retval$stderr[1])
            }
            retval$stdout <- suppressWarnings(fread(stdoutFile))
            retval$stderr <- readLines(stderrFile)
            unlink(c(stdoutFile, stderrFile))
            retval$retry <- FALSE
            return(retval)
        }, error =  function(e) {
            if (restart_times > 0){
                return(list(retry = TRUE))
                message("Call failed:", e, " retrying(",restart_times," tries left)")                
            } else {
                stop(e)
            }
            Sys.sleep(retry_wait)
        })
        retry <- tmp$retry
        restart_times  <- restart_times - 1
    }
    return(tmp)
}

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
    status <- "@should_never_be returned"
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
            s <- "@PASS"
        } else {
            s <- "@insufficient_neutral_mutations"
        }        
    }
    return(list(out = out, s = s))
}

## load global model
load(args$global)
global_model <- model



## Load in blocks
blocks <- fread(args$blocks, header = FALSE)

pass_count <- 0
for(i in seq_len(nrow(blocks))) {
    # If just starting iterating over bed file overwrite existing output files
    if(i == 1) {
        append_out <- FALSE
    } else {
        append_out <- TRUE
    }
    
    chr <- blocks[i, ]$V1
    bin_start <- blocks[i, ]$V2
    bin_end <- blocks[i, ]$V3
   
    # Add flanking regions to co-ordinates
    start <- max(0, bin_start - args$flanking_size) # if start is less than 0 make it 0
    end <- bin_end + args$flanking_size
    s <- "@should_not_be_seen"
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
    str = paste("tabix {mutation_file} {chr}:{start}-{end} | perl filter_mutation.pl {min_coverage}",
        "{max_frequency} | perl add_context_mutation_rate.pl {context_frequency_table} |",
        "bedtools intersect -a - -b {covariate_file} -sorted -wo | cut -f1-7,11- |",
        "intersectBed -loj -a - -b {neutral_file} -sorted")
    cmd <- strformat(str, args_local)
    ## Run command and read into data.table
    d <- suppressWarnings(robust_system(cmd, 10, 5))
    ## Check that data read in as expected
    if(d$exitStatus != 0){
        s <- paste0("@variant_retrieval_exit_status(",d$exitStatus,")")
        stop(d$stderr, "(", d$exitStatus,")")
    } else if(nrow(d$stdout) == 0) {
        s <- "@no_regional_post_filter_mutations"
    } else if(ncol(d$stdout) != 13) {
        s <- paste0("@malformed_input_ncol_", ncol(d))
    } else {
        ## Rename columns
        setnames(d$stdout, colnames(d$stdout),
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
        results <- compute_local_mutation_rates(d$stdout, args$min_n_mutation, global_model)
        s <- results$s
    }
    if (s == "@PASS") {
        ## Remove from the output all regions outside the original query block, remember we added in some
        ## flanks to increase the number of sites to compute the local adjustment values
        results$out <- results$out[start >= bin_start & end < bin_end]
        fwrite(results$out, file = args$output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
               append = append_out, compress = "auto")
        pass_count <- pass_count + 1
    }
    blocks[i, status := s]
    fwrite(blocks[i, ], file = args$log, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
           compress = "auto", append = append_out)
    ## A bit of cleanup
    rm(results)
}

## Always write end of file line to deal with snakemake complaint of missing file if no blocks pass
fwrite(data.table("#END_OF_FILE"), file = args$output, quote=FALSE, sep="\t", row.names=FALSE,
       col.names=FALSE, append = TRUE, compress = "auto")

