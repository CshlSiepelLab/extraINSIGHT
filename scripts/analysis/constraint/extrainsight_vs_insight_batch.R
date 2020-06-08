#! /usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("-b","--bed-dir", dest="bed_dir", type = "character",
                    required=TRUE, help="directory containing bed files of query regions")
parser$add_argument("-g","--bed_genome", dest="bed_genome", type="character",
                    required=TRUE,
                    choices=c('hg19', 'hg38'),
                    help="genome version the bed files are from. Must all be same. Either hg19 or hg38")
parser$add_argument("--insight-genome", dest="insight_genome", type="character", default="hg19",
                    choices=c('hg19'),
                    help="The genome version that INSIGHT is using. Do not touch unless you are sure.")
parser$add_argument("--extrainsight-genome", dest="extra_insight_genome", type="character", default="hg38",
                    choices=c('hg19', 'hg38'),
                    help="The genome version that ExtRaINSIGHT is using. Do not touch unless you are sure.")
parser$add_argument("-o","--out-dir", dest="out_dir", type="character",
                    default = ".", help="results output directory (default = '.')")
parser$add_argument("--exclude-chr", dest="exclude_chr", nargs="+", type="character",
                    required=FALSE, help="chromosomes to exclude from the analysis")
parser$add_argument("--max-sites", dest="max_sites", type="integer", default = 1e5 ,
                    required=FALSE, help="maximum number of sites to analyze, if exceeded it will subsample")
parser$add_argument("--max-processes", dest="max_processes", type="integer", default = 1 ,
                    required=FALSE, help="maximum number of local jobs to run at once")
args <- parser$parse_args()

library(processx)

## tst <- "-b ../../../data/grch38/annotation_bed/top_5_expressed -o ../../../results/grch38/gnomad_v3.0/constraint/top_5_expressed --max-sites 10000 -g hg38 --extrainsight-genome hg38 --max-processes 5 --exclude-chr X"
## args <- parser$parse_args(unlist(strsplit(tst, "\\s+")))

args$bed_dir <- normalizePath(args$bed_dir)
args$out_dir <- normalizePath(args$out_dir)
args$log_dir <- file.path(args$out_dir, "log")

dir.create(args$out_dir, FALSE, TRUE)
dir.create(args$log_dir, FALSE, TRUE)

## Gets the location of this_file
## Lifted from https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
this_file <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

## Mimics pythons f-strings
## Function lifted from https://stackoverflow.com/questions/44763056/is-there-an-r-equivalent-of-pythons-string-format-function
strformat <- function(str, vals) {
    vars <- stringi::stri_match_all(str, regex = "\\{(.*?)\\}", vectorize_all = FALSE)[[1]][,2]
    x <- str
    for (i in seq_along(vars)) {
        varName <- vars[i]
        varCode <- paste0("{", varName, "}")
        if(!is.null(vals[[varName]])) {
            x <- stringi::stri_replace_all_fixed(x, varCode, vals[[varName]], vectorize_all = TRUE)
        } else {
            warning("Missing variable in vals: ", varName )
        }
    }
    return(x)
}

## Count number of live processx processes given a list of processes
live_processes <- function(z) {
    return(sum(unlist(lapply(z, function(x) x$is_alive()))))
}

## Get bed files in bed directory
bed_files <- dir(args$bed_dir)[grepl(".bed", dir(args$bed_dir))]

## Creates template string for optional exclude chromosome argument
exclude_string <- ""
if (!is.null(args$exclude_chr)) {
    exclude_string <- "--exclude-chr {exclude_chr}"
}

## Setup command that specific instances will be generated of
script <- file.path(dirname(this_file()), "insight_vs_extraINSIGHT.py")
cmd_args_template <- paste("-b {bed_dir}/{bed} -g {bed_genome} --extrainsight-genome",
                           "{extra_insight_genome} -o {out_dir}/{core}",
                           "--insight-genome {insight_genome}",
                           "--max-sites {max_sites}", exclude_string)
out_template <- "{log_dir}/{core}.out"
err_template <- "{log_dir}/{core}.err"

proc_list <- list()
## Iterate over bed files and dispatch one instance of a command per bed file
for(b_ind in seq_along(bed_files)){
    print(paste("Dispatching job:", bed_files[b_ind]))
    b <- bed_files[b_ind]
    core_name <- gsub("\\.bed.*","",b)
    ## Fill in template strings with specific configurations
    tmp_list <- c(args, list(bed = b, core = tolower(core_name)))
    cmd_inst <- strformat(cmd_args_template, tmp_list)
    out_inst <- strformat(out_template, tmp_list)
    err_inst <- strformat(err_template, tmp_list)
    ## Launch process
    p <- process$new(script, args = unlist(strsplit(cmd_inst, "\\s+")),
                     stdout = out_inst, stderr = err_inst)
    proc_list <- append(proc_list, p)
    ## Monitor number of processes and launch new one when there are fewer
    ## than the max permissible number of processes
    while (live_processes(proc_list) >= args$max_processes) {
        Sys.sleep(10)
    }
}


