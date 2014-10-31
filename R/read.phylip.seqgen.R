read.phylip.seqgen <- function(text, byrow = TRUE, code.type = "NUCLEOTIDE") {
    ######## code adapted from read.seqgen in package 'phyclust' by Chen and Dorman
    phylip <- list(code.type = "NUCLEOTIDE", info = NULL, nseq = NULL, seqlen = NULL, 
        seqname = NULL, org.code = NULL, byrow = byrow)
    
    phylip$info <- text[1]
    
    splitter <- function(x, y, ...) {
        unlist(strsplit(x, y, ...))
    }
    
    tmp <- splitter(phylip$info, " ")
    tmp <- tmp[tmp != ""]
    phylip$nseq <- as.numeric(tmp[1])
    phylip$seqlen <- as.numeric(tmp[2])
    
    #### splitting names from sequences this way will only work if there are never
    #### spaces in names.
    nameseq <- do.call("rbind", lapply(text[2:length(text)], splitter, " "))
    org.names <- nameseq[, 1]
    org.names <- gsub("  *", "", org.names)
    org <- do.call("rbind", lapply(nameseq[, 2], splitter, ""))
    org <- as.matrix(org, nrow = phylip$nseq)
    
    phylip$seqname <- org.names
    if (code.type == "NUCLEOTIDE") {
        phylip$org.code <- matrix(org, nrow = phylip$nseq, ncol = phylip$seqlen)
    } else {
        stop("The code.type is not found.")
    }
    if (!byrow) {
        phylip$org.code <- t(phylip$org.code)
    }
    phylip$byrow <- byrow
    class(phylip) <- "seq.data"
    phylip
} 
