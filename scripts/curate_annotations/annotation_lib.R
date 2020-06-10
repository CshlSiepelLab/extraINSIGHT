library(GenomicRanges)
library(rtracklayer)
library(stringi)

export_bed <- function(gr, file_name, path = ".", adjust_start = 0, adjust_end = 0, style = "UCSC",
                       sort_bed= TRUE) {
    dir.create(path, FALSE, TRUE)
    start(gr) <- start(gr) + adjust_start
    end(gr) <- end(gr) + adjust_end
    out_file <- file.path(path, file_name)
    GenomeInfoDb::seqlevelsStyle(gr) <- style
    if (sort_bed) {
        seqlevels(gr) <- sort(as.character(seqlevels(gr)))
        gr <- sort(gr)
    }
    if(tools::file_ext(file_name) == "gz") {
        con <- gzfile(out_file)
    } else {
        con <- gzfile(out_file)
    }
    rtracklayer::export.bed(gr, con)
}

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
