library(ensembldb)
library(AnnotationHub)
## BiocManager::install("EnsDb.Hsapiens.v75")

## This function will return the appropriate txdb
## hg19 is from a pre-built package
## grch38 is from a pre-built annotation hub
load_txdb <- function(id) {
    if (id == "GRCh38") {
        hub <- AnnotationHub()
        txdb <- query(hub, c("EnsDb", "sapiens", "100"))[[1]]
    } else if (id == "GRCh37") {
        library(EnsDb.Hsapiens.v75)
        txdb <- EnsDb.Hsapiens.v75
    } else {
        error("unsupported id")
    }
    return(txdb)    
}
