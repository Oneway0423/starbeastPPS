gene.tree.prob <- function(sptree, gtree, association, ploidy = 1) {
    
    species <- sptree$tip.label
    nspec <- length(species)
    nbranch <- (2 * nspec) - 1
    
    tiplist <- list()
    for (i in 1:nspec) {
        tiplist[[species[i]]] <- association[which(association[, 1] == species[i]), 
            2]
    }
    
    #### this subfunction returns the tip numbers for a given node
    node.tips <- function(phy, node) {
        n <- length(phy$tip.label)
        if (node <= n) {
            node
        } else {
            l <- numeric()
            d <- phy$edge[which(phy$edge[, 1] == node), 2]
            for (j in d) {
                if (j <= n) {
                  l <- c(l, j)
                } else {
                  l <- c(l, node.tips(phy, j))
                }
            }
            l
        }
    }
    ##### end subfunction subfunction that gets a tree containing descendants of all gene
    ##### lineages that pass through a branch do I really need that node conditional???
    ##### or should I just deal with it later?
    get.lin <- function(sptree, gtree, association, node) {
        species <- sptree$tip.label
        if ((length(species) + 1) != node) {
            spec.desc <- species[node.tips(sptree, node)]
            gtree.desc <- c()
            for (i in 1:length(spec.desc)) {
                gtree.desc <- c(gtree.desc, association[which(association[, 1] == 
                  spec.desc[i]), 2])
            }
            index <- match(gtree.desc, gtree$tip.label)
            pruned <- drop.tip(gtree, gtree$tip.label[-index])
        } else {
            pruned <- gtree
        }
        return(pruned)
    }
    ##### end subfunction prepare branch demographic info...
    demo.maker <- function(sptree, nspec) {
        demo <- cbind(sptree$edge[, 2], sptree$edge.length)
        demo <- rbind(demo, c((nspec + 1), Inf))
        demo <- cbind(demo, sptree$dmv)
        demo <- demo[order(demo[, 1]), ]
        sbt <- branching.times(sptree)
        sbt <- sbt[order(as.numeric(names(sbt)))]
        sbt <- c(rep(0, nspec), sbt)
        demo <- cbind(demo, sbt)
        colnames(demo) <- c("node", "length", "dmv", "sbt")
        rownames(demo) <- c(1:length(sbt))
        return(demo)
    }
    ##### 
    
    ##### now how to calculate likelihoods!??!?! for a single branch...
    
    b.prob <- function(sptree, gtree, demo, node) {
        pruned <- get.lin(sptree, gtree, association, node)
        gbt <- sort(branching.times(pruned))
        gbt <- c(0, gbt)
        gbt <- cbind(gbt, length(gbt):1)
        start <- demo[node, "sbt"]
        end <- (demo[node, "sbt"] + demo[node, "length"])
        enter <- gbt[gbt[, 1] == max(gbt[gbt[, 1] <= start, 1]), 2][1]  ## need to double check that adding [1]
        exit <- gbt[gbt[, 1] == max(gbt[gbt[, 1] <= end, 1]), 2][1]  ##   actually does what's expected!
        
        #### this 'if' is made irrelevant by ignoring one-tip species below.
        if (start == 1) {
            return(1)
        } else {
            gbt <- gbt[gbt[, 1] >= start & gbt[, 1] <= end, ]
            
            if (start != 0) {
                gbt <- rbind(c(start, enter), gbt, c(end, exit))
            } else {
                gbt <- rbind(gbt, c(end, exit))
            }
            
            waits <- gbt[2:length(gbt[, 1]), 1] - gbt[1:length(gbt[, 1]) - 1, 1]
            waits <- cbind(waits, gbt[1:length(gbt[, 1]) - 1, 2])
            lambda <- ((waits[, 2] * (waits[, 2] - 1))/((2 * demo[node, "dmv"]) * 
                ploidy))
            exponent <- exp(-lambda * waits[, 1])
            
            exponent <- prod(exponent[lambda != 0])
            lambda <- prod(lambda[1:(length(lambda) - 1)])
            probability <- lambda * exponent
        }
        
        
        return(probability)
    }
    
    ###### now actually calculate the whole gene tree probability remember to exclude
    ###### species tip branches with only one allele sampled
    
    demo <- demo.maker(sptree, nspec)
    numtips <- lapply(tiplist, length)
    if (any(numtips == 1)) {
        demo2 <- demo[-which(numtips == 1), ]
    } else {
        demo2 <- demo
    }
    
    lnP <- c()
    nodes <- demo2[, "node"]
    names(nodes) <- NULL
    for (i in 1:length(nodes)) {
        lnP <- c(lnP, log(b.prob(sptree, gtree, demo, nodes[i])))
    }
    
    return((sum(lnP)))
} 
