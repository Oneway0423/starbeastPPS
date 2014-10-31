analyze.coalescent <- function(data, msdir = "/directory/containing/ms") {
    ntrees <- length(data$genes)
    nsamples <- length(data$species.trees)
    
    ### data structures
    
    empirical.gene.coal <- matrix(nrow = nsamples, ncol = ntrees)
    empirical.gene.probs <- matrix(nrow = nsamples, ncol = ntrees)
    colnames(empirical.gene.probs) <- data$genes
    colnames(empirical.gene.coal) <- data$genes
    ## empirical.sum.probs<-c()
    
    
    ## simulated.genes<-list()###should I actually bother to save the gene trees?
    simulated.gene.coal <- matrix(nrow = nsamples, ncol = ntrees)
    simulated.gene.probs <- matrix(nrow = nsamples, ncol = ntrees)
    colnames(simulated.gene.probs) <- data$genes
    colnames(simulated.gene.coal) <- data$genes
    ## simulated.sum.probs<-c()
    
    #### calculate empirical gene tree stats
    for (j in 1:ntrees) {
        for (i in 1:nsamples) {
            empirical.gene.probs[i, j] <- gene.tree.prob(data$species.trees[[i]], 
                data$gene.trees[[j]][[i]], data$associations[["empirical"]][[j]], 
                ploidy = data$ploidy[j])
            empirical.gene.coal[i, j] <- deep.coal(data$species.trees[[i]], data$gene.trees[[j]][[i]], 
                data$associations[["empirical"]][[j]])
        }
        cat("empirical tree ", j, " is done!\n")
    }
    #### simulate new gene trees and calculate their stats
    for (j in 1:ntrees) {
        ns <- length(data$gene.trees[[j]][[1]]$tip.label)
        for (i in 1:nsamples) {
            simtree <- mstree(data$species.trees[[i]], msdir, nseq = ns, nreps = 1, 
                samplescheme = data$associations[["simulate"]][[j]], ploidy = data$ploidy[j])
            simulated.gene.probs[i, j] <- gene.tree.prob(data$species.trees[[i]], 
                simtree, data$associations[["calculate"]][[j]], ploidy = data$ploidy[j])
            simulated.gene.coal[i, j] <- deep.coal(data$species.trees[[i]], simtree, 
                data$associations[["calculate"]][[j]])
        }
        cat("simulated tree ", j, " is done!\n")
    }
    
    result.all <- list()
    result.prob <- list()
    result.coal <- list()
    result.prob[["empirical"]] <- empirical.gene.probs
    result.prob[["simulated"]] <- simulated.gene.probs
    result.prob[["test.stat"]] <- simulated.gene.probs - empirical.gene.probs
    result.coal[["empirical"]] <- empirical.gene.coal
    result.coal[["simulated"]] <- simulated.gene.coal
    result.coal[["test.stat"]] <- simulated.gene.coal - empirical.gene.coal
    
    result.all[["probs"]] <- result.prob
    result.all[["coal"]] <- result.coal
    class(result.all) <- "coalescentteststats"
    return(result.all)
} 
