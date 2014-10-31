read.gtree.samples <- function(file) {
    ##### parts of this function were adapted from a function in package 'phyloch',
    ##### written by Christoph Heibl read in trees without stats

    ##### get tree strings
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    X <- X[grep("tree STATE", X)]
    X <- gsub("tree STATE_.*\\[&R\\] ", "", X)

    ###### a function to produce branch rate vectors and order them appropriately

    extract.order.stats <- function(treestring) {

        ntax <- (length(unlist(strsplit(treestring, "\\["))) + 1)/2
        edges <- (ntax + 2):((2 * ntax) - 1)

        Y <- treestring


        #### adds internal node labels
        for (i in 1:length(edges)) {
            repl <- paste(")", edges[i], ":[", sep = "")
            Y <- sub("):\\[", repl, Y)
        }

        meta <- unlist(strsplit(Y, "\\[|\\]"))[grep("rate", unlist(strsplit(Y, "\\[|\\]")))]

        metacols <- length(unlist(strsplit(meta[1], ",")))
        meta <- gsub("&|rate=|\\{|\\}", "", meta)

        Z <- gsub(";", "", Y)
        Ysub <- gsub("\\[[^]]*\\]", "\\[\\]", Z)
        Ysub <- unlist(strsplit(Ysub, ",|)"))
        Ysub <- gsub("\\(|\\)|;|\\[|\\]", "", Ysub)

        brlen <- do.call("rbind", strsplit(Ysub, ":"))
        brrate <- do.call("rbind", strsplit(meta, ","))
        branchdata <- cbind(brlen, brrate)

        if (metacols == 1) {
            colnames(branchdata) <- c("br", "length", "rate")
        }
        rownames(branchdata) <- branchdata[, 1]

        string <- gsub("\\[[^]]*\\]", "", Y)

        stree <- read.tree(text = string)
        translate <- cbind(stree$node.label[-1], (ntax + 2):edges[length(edges)])
        translate <- rbind(translate, cbind(stree$tip.label, 1:ntax))

        rownames(translate) <- translate[, 2]
        translate2 <- translate[as.character(stree$edge[, 2]), ]
        branchdata2 <- branchdata[translate2[, 1], ]

        rownames(branchdata2) <- NULL
        list(tree=stree, branchdata=branchdata2)

    }
    tree_info <- lapply(X, extract.order.stats)
    tree <- lapply(tree_info, function(x) x$tree)
    branchdata <- lapply(tree_info, function(x) x$branchdata)

    ####### append branch rate stats to trees.

    if (length(tree) > 1) {
        tree <- mapply(function(tr, mat) {
            mat <- mat[, c(-1, -2), drop=FALSE]
            mat <- apply(mat, 2, function(x) {
                as.numeric(gsub("(.+)\\.([0-9]+\\.[0-9]+E?-?[0-9]?)$",
                                "\\2", x))
            })
            tr[["rate"]] <- mat
            tr
        }, tr=tree, mat=branchdata, SIMPLIFY=FALSE)
    } else {
        mode(branchdata[[1]]) <- "numeric"
        tree[["rate"]] <- as.numeric(branchdata[[1]][, c(-1, -2)])
    }
    return(tree)
}
