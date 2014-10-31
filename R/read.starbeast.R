read.starbeast <- function(beast.xml, combinedfiledirectory, logfile = "combined.log") {
    
    ##### read in a beast xml file in order to extract the taxon-allele bindings and
    ##### locus names.
    X <- scan(file = beast.xml, what = "", sep = "\n", quiet = TRUE)
    X <- gsub("\t", "", X)
    
    Y <- X[grep("<alignment id=\".*\" dataType=\"nucleotide\">", X)[1]:grep("</alignment>", 
        X)[length(grep("</alignment>", X))]]
    Q <- X[grep("<treeLikelihood id=", X)[1]:grep("</treeLikelihood>", X)[length(grep("</treeLikelihood>", 
        X))]]
    
    X <- X[(grep("<species id=\"species\">", X) + 1):(grep("</species>", X) - 1)]
    
    species <- X[grep("<sp id=", X)]
    species <- gsub("<sp id=\"", "", species)
    species <- gsub("\">", "", species)
    
    alleles <- X[grep("<taxon idref=", X)]
    alleles <- gsub("<taxon idref=\"", "", alleles)
    alleles <- gsub("\"/>", "", alleles)
    
    num <- (grep("</sp>", X) - 1) - grep("<sp id=", X)
    
    spec <- c()
    for (i in 1:length(species)) {
        spec <- c(spec, rep(species[i], times = num[i]))
    }
    
    bindings <- cbind(spec, alleles)
    
    loci <- X[(grep("<!-- Collection of Gene Trees", X) + 2):(grep("</geneTrees>", 
        X) - 1)]
    
    if (any(grep("ploidy", loci))) {
        ploidy <- loci[grep("ploidy", loci)]
        ploidy <- gsub("<gtree ploidy=\"", "", ploidy)
        ploidy <- gsub("\">", "", ploidy)
    } else {
        ploidy <- rep(1, length(grep("<treeModel idref=", loci)))
    }
    
    loci <- X[grep("<treeModel idref=", X)]
    loci <- gsub("<treeModel idref=\"", "", loci)
    loci <- gsub(".treeModel\"/>", "", loci)
    
    ### reading sequence alignments...
    
    p.begin <- grep("<treeLikelihood id=", Q)
    p.end <- grep("</treeLikelihood>", Q)
    genes.to.partitions <- matrix(ncol = 2, nrow = length(p.begin))
    colnames(genes.to.partitions) <- c("genes", "partitions")
    for (i in 1:length(p.begin)) {
        element <- Q[p.begin[i]:p.end[i]]
        gene <- element[grep("treeModel", element)]
        gene <- gsub("<treeModel idref=\"", "", gene)
        gene <- gsub(".treeModel\"/>", "", gene)
        part <- element[grep("pattern", element)]
        part <- gsub("<patterns idref=\"", "", part)
        part <- gsub(".patterns\"/>", "", part)
        genes.to.partitions[i, ] <- c(gene, part)
    }
    
    alignments <- list()
    a.begin <- grep("<alignment id=\".*\" dataType=\"nucleotide\">", Y)
    a.end <- grep("</alignment>", Y)
    
    for (i in 1:length(genes.to.partitions[, 1])) {
        Z <- Y[a.begin[i]:a.end[i]]
        
        ID <- Z[grep("<taxon idref=", Z)]
        ID <- gsub("<taxon idref=\"", "", ID)
        ID <- gsub("\"/>", "", ID)
        
        seq <- Z[(grep("<taxon idref=", Z) + 1)]
        # seq <- as.list(seq)
        names(seq) <- ID
        
        seq <- sapply(seq, strsplit, split = "")
        
        alignments[[genes.to.partitions[i, 2]]] <- seq
        
    }
    
    cat("bindings, loci, and alignments read\n")
    
    ### read in gene trees from combinedfiledirectory. assumes custom perl script used
    ### to thin and burn in beast results.
    gene.trees <- list()
    for (i in 1:length(loci)) {
        treefile <- c(combinedfiledirectory, "/", loci[i], ".combined.trees")
        gene.trees[[loci[i]]] <- read.gtree.samples(paste(unlist(treefile), collapse = ""))
        cat("locus ", loci[i], " is done.\n")
    }
    cat("gene trees done\n")
    
    ### read in species trees from combinedfiledirectory. assumes custom perl script
    ### used to thin and burn in beast results.
    species.trees <- read.sptree.samples(paste(unlist(c(combinedfiledirectory, "/", 
        "species.combined.trees")), collapse = ""))
    cat("species trees done\n")
    
    ### make the various taxon associations.
    empirical <- list()
    for (i in 1:length(loci)) {
        locus.binding <- c()
        tips <- gene.trees[[loci[[i]]]][[1]]$tip.label
        for (j in 1:length(tips)) {
            locus.binding <- rbind(locus.binding, bindings[which(bindings[, 2] == 
                tips[j]), ])
        }
        empirical[[loci[[i]]]] <- as.matrix(as.data.frame(locus.binding))
    }
    
    simulate <- list()
    for (i in 1:length(loci)) {
        simulate[[loci[[i]]]] <- as.matrix(as.data.frame(table(empirical[[loci[[i]]]][, 
            1])))
        # simulate[[loci[[i]]]][,1]<-levels(simulate[[loci[[i]]]][,1])
    }
    
    calculate <- list()
    for (i in 1:length(loci)) {
        spec <- c()
        for (j in 1:length(simulate[[i]][, 1])) {
            spec <- c(spec, rep(simulate[[i]][j, 1], times = simulate[[i]][j, 2]))
        }
        calculate[[loci[[i]]]] <- cbind(spec, 1:length(spec))
        calculate[[loci[[i]]]] <- as.matrix(as.data.frame(calculate[[loci[[i]]]]))
    }
    
    associations <- list()
    associations[["empirical"]] <- empirical
    associations[["simulate"]] <- simulate
    associations[["calculate"]] <- calculate
    
    #### read in log file
    logfile <- c(combinedfiledirectory, "/", logfile)
    logdata <- read.table((paste(unlist(logfile), collapse = "")), header = TRUE)
    cat("logfile done\n")
    #### organize output components
    output <- list()
    output[["associations"]] <- associations
    output[["alignments"]] <- alignments
    output[["gene.trees"]] <- gene.trees
    output[["species.trees"]] <- species.trees
    output[["log"]] <- logdata
    output[["genes"]] <- loci
    output[["genes.to.partitions"]] <- genes.to.partitions
    output[["ploidy"]] <- as.numeric(ploidy)
    class(output) <- "starbeast.data"
    return(output)
    
} 
