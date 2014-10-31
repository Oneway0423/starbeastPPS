get.model <- function(data, partition) {
    
    locpar <- names(data$log)[grep(paste(partition, "\\.", sep = ""), names(data$log))]
    
    ## base frequencies
    if (any(grep("frequencies", locpar))) {
        basefreq = "estimated"
    } else {
        basefreq = "fixed"
    }
    
    ## rate matrix
    if (any(grep("kappa1", locpar))) {
        rmat = "trn"
    } else {
        if (any(grep("kappa", locpar))) {
            rmat = "hky"
        } else {
            rmat = "gtr"
        }
    }
    
    ## among-site rate variation
    if (any(grep("pInv", locpar))) {
        pInv <- TRUE
    } else {
        pInv <- FALSE
    }
    if (any(grep("alpha", locpar))) {
        alpha <- TRUE
    } else {
        alpha <- FALSE
    }
    
    model <- list()
    model[["rmat"]] <- rmat
    model[["basefreq"]] <- basefreq
    model[["pInv"]] <- pInv
    model[["alpha"]] <- alpha
    
    return(model)
} 
