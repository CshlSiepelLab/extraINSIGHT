library(data.table)

read_insight_model_block <- function(block_header, file) {
    start <- grep(block_header,readLines(file))
    section_ends <- grep("^$",readLines(file))
    end <- min(section_ends[section_ends > start])
    out <- fread(file, skip = start, nrows = end - start - 1)
    out[, V1 := tolower(gsub(":","",V1))]
    setnames(out, c("V1","V2"), c("name", "val"))
    return(out)
}

extract_conservation_estimates <- function(directory_list, residual_rho = FALSE){
    out <- list()
    for(d in directory_list) {
        insight <- file.path(d, "INSIGHT","i_hg19.model")
        ei <- file.path(d, "ExtRaINSIGHT","strong_selection_estimate.txt")
        if(!file.exists(insight) || !file.exists(ei)){
            warning("Missing one or more output files in:", d)
            next
        }
        i_param <- read_insight_model_block("PARAMETERS", insight)[name == "rho",]
        i_se <- read_insight_model_block("UNCERTAINTY", insight)[name == "rho",]
        i_param$se <- i_se$val
        ss <- fread(ei)
        ss_param <- ss[name == "strong_selection"]
        ss_se <- ss[name == "se"]
        ss_param$se <- ss_se$val
        if (residual_rho) {
            i_param$val<- i_param$val - ss_param$val
        }
        out[[basename(d)]] <- rbind(i_param, ss_param)
    }
    final <- rbindlist(out, idcol = "label")
    setnames(final, "val", c("absolute", "residual")[residual_rho + 1])
    return(final)
}
