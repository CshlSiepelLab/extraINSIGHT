library(argparse)

parser <- ArgumentParser()
parser$add_argument("-m","--mutation-file", dest="mutation_file", type="character",
                    required=TRUE, help="file containing per allele sitewise mutation probabilities")
parser$add_argument("-n","--neutral", dest="neutral_file", type="character",
                    required=TRUE, help="neutral regions in bed format")
parser$add_argument("--num-sites", dest="num_sites", type="integer",
                    required=FALSE, help="number of sites to subsample for analysis")
parser$add_argument("-o", "--out-dir", dest="out_dir", type="character",
                    required=FALSE, help="directory to write output plots to", default = ".")
parser$add_argument("--fraction-sites", dest="percent_sites", type="double",
                    required=FALSE, help="fraction of sites to subsample for analysis")


args <- parser$parse_args()
## tst <- "-m ../../../results/grch38/gnomad_v3.0/mutation_model/final_mutation_rates.bed.gz -n ../../../data/grch37/neutral_regions/hg19_neutral_region.bed --num-sites 100000"
## tst <- unlist(strsplit(tst, " "))
## args <- parser$parse_args(tst)

library(IRanges)
library(data.table)
library(ggplot2)
library(cowplot)

# Read-in neutral regions
bed <- fread(args$neutral)

# If both number and percent of sites null, use all sites
if (is.null(args$num_sites) &&  is.null(args$percent_sites)) {
    args$num_sites <- sum(unlist(bed[,3] - bed[, 2]))
} else if (!is.null(args$percent_sites)) {
    args$num_sites <- round(sum(unlist(bed[,3] - bed[, 2])) * args$percent_sites)
}

# Function lifted from https://stackoverflow.com/questions/44763056/is-there-an-r-equivalent-of-pythons-string-format-function
strformat <- function(str, vals) {
    vars <- stringi::stri_match_all(str, regex = "\\{(.*?)\\}", vectorize_all = FALSE)[[1]][,2]
    x <- str
    for (i in seq_along(vars)) {
        varName <- vars[i]
        varCode <- paste0("{", varName, "}")
        x <- stringi::stri_replace_all_fixed(x, varCode, vals[[varName]], vectorize_all = TRUE)
    }
    return(x)
}

subsample_bed <- function(bed, n, sort_bed = FALSE) {
    setnames(bed, colnames(bed)[1:3], c("chrom", "start", "end"))
    if (sort_bed) {
        bed <- bed[order(chrom, start, end)]
    }
    region_lengths <- bed$end - bed$start
    total_length <- sum(region_lengths)
    ## Choose which sites to keep
    sites <- sort(sample(seq_len(total_length), size = n, replace = FALSE))
    ## Compute cumulative length distribution
    cum_lengths <- cumsum(region_lengths)
    ## Identify sites for each region
    groups <- S4Vectors::subjectHits(
        IRanges::findOverlaps(sites,
                              IRanges::IRanges(cum_lengths - region_lengths + 1, cum_lengths))
        )
    chr <- bed[groups, ]$chrom
    start <- sites - cum_lengths[groups] + region_lengths[groups] + bed[groups, ]$end -1
    out <- data.table(chrom = chr, start = start, end = start + 1)
    setkeyv(out, colnames(out))
    return(out)
}

## Subsample sites and write to tempfile
sub_sites <- subsample_bed(bed, args$num_sites)
tab <- copy(sub_sites)
tab[,start:=start+1]
tmp_file_bed <- paste0(tempfile(),".bed.gz")
tmp_file_tab <- paste0(tempfile(),".tab.gz")
fwrite(sub_sites, tmp_file_bed, col.names = FALSE, sep = "\t")
fwrite(tab, tmp_file_tab, col.names = FALSE, sep = "\t")

## Create tabix query
cmd <- strformat(str = "tabix {mutation_file} -R {tmp_file}",
                 vals = c(args,list(tmp_file = tmp_file_bed)))
bed_in <- fread(cmd = cmd)
setnames(bed_in, colnames(bed_in), c("chrom",
                                     "start",
                                     "end",
                                     "context",
                                     "mutation_flag",
                                     "coverage",
                                     "mutation_rate")
)

## Cut mutation probabilities into bins
brks <- seq(0,1,0.001)
bed_in[, bin := findInterval(mutation_rate, brks, all.inside = T)]
ove <- bed_in[,mean(mutation_rate), by = "bin"]
ove[, midpoint:= (brks[bin] + brks[bin + 1]) / 2]

g <- ggplot(ove, aes(x= midpoint, y=V1)) + 
    geom_point()+
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red")+
    theme_cowplot()+
    xlab("Probability of mutation")+
    ylab("Observed fraction of mutated sites")

dir.create(args$out_dir,FALSE,TRUE)
ggsave(filename = "mutation_model_calibration.pdf", path=args$out_dir, plot = g)
