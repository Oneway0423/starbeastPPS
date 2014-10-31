gap.columns <- function(gapmatrix) {
    dims <- dim(gapmatrix)
    nonmissing.rows <- rowSums(gapmatrix)
    nonmissing.rows <- nonmissing.rows < dims[2]
    gap.mat <- gapmatrix[nonmissing.rows, ]
    gap.cols <- colSums(gap.mat) == 0
    return(gap.cols)
} 
