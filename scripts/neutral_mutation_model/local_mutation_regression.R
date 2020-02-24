library(data.table)
args = commandArgs(trailingOnly=TRUE)

## load global model
load(args[1])
 
min_n_mutation = as.numeric(args[2])

## Read in local data from stdin
d <- fread('cat /dev/stdin')

## Check that data read in
if(!exists("d")){
    stop("Failed to read in mutation data")
}

if(ncol(d) != 0){
    message("No neutral mutations passed the filters for this region, exiting...")
    # write("NULL", file = stdout())
    q("no")
}

if(ncol(d) != 12){
    message("Malformed input (top row printed):")
    message(head(d, 1))
    stop(paste("Something has changed in data processing, input has",ncol(d),"columns, should have 12. Has a covariate been added?"))
}

## It should have 12 columns which are:
## 1)chrom
## 2)start
## 3)end
## 4)context-alternate_allele
## 5)mutation_flag (0 = no mutation, 1 = mutation at that particular site)
## 6)average sequencing coverage
## 7)log(mutation probability) based on context from mutation_frequency_table.txt
## 8)GC_content
## 9)number of bases that this site overlaps with the the covariate bed file by
## 10)chrom (of overlapping neutral region, or "." if not in neutral region) 
## 11)start (of overlapping neutral region, or "-1" if not in neutral region)
## 12)end (of overlapping neutral region, or "-1" if not in neutral region)

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
           "base_overlap",
           "chrom_neutral",
           "start_neutral",
           "end_neutral"
           )
         )

out_cols = c(colnames(d)[1:6], "mutation_scaling_factor")

## Now check that there are enough sites for which there are mutations and are not overlapped with neutral regions
mutation_count <- sum(d$mutation_flag== 1 & d$chrom_neutral!= '.')
if (mutation_count >= min_n_mutation){
    ## Get the predicted mutation rate for each site based on the global model
    d[, predicted_mutation_rate := predict(model, newdata=d)]
    ## Subset out the neutral sites
    d_neut <- d[chrom_neutral != ".", ]
    ## Compute the scaling factor necessary to make the global model predictions fit the observed local neutral
    ## mutation rate as well as possible    
    local_model <- glm(mutation_flag~predicted_mutation_rate, data=d_neut, family="binomial")
    d[,mutation_scaling_factor:= predict(local_model, newdata=d, type="response")]
    ## Write out the data
    write.table(d[, ..out_cols], stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
} else {
    message(paste0("Insufficient number of mutations passed the filters(", mutation_count,") ",min_n_mutation," required"))
}

