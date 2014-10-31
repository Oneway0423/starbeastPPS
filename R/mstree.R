##' takes a species tree of object 'phylo' with a vector giving branch
##' effective population sizes and simulates a coalescent
##' genealogy. this is used internally, but may be useful for some
##' people so it is documented here.
##'
##' simulates coalescent genealogies from species trees. effective
##' population sizes on each branch must be specified in a 'dmv'
##' vector appended to the 'phylo' object.
##'
##' ms should be to increase the number of digits in the branch
##' calculations.
##'
##' @title Simulate Coalescent Genealogies Given a Species Tree using Hudson's ms
##' @param phy the population tree in which gene trees are to be
##' simulated. IN THIS SCRIPT, ALL TREES MUST BE SPECIES TREES WITH
##' ASSOCIATED 'DMV' BRANCH WIDTH VALUES
##' @param msdir the directory containing Hudson's ms.
##' @param nseq the number of alleles (tips in the coalesecent
##' genealogies) to simulate
##' @param nreps the number of trees to simulate
##' @param samplescheme a data frame with 2 columns and nseq rows. the
##' first column contains the tip labels of the species tree and the
##' second contains the number of alleles to be sampled from the
##' corresponding population.
##' @param ploidy a multiplier for the 'dmv' vector. if ploidy is
##' specified, then dmv=dmv*ploidy. do not specify it unless you wish
##' the dmv vector to be scaled by ploidy. refers to *BEAST's 'ploidy'
##' in the xml files and modifies DMV values. When all loci have the
##' same ploidy, it should left as 1. When ploidy varies, it should be
##' 0.5 mitochondrial and 2 for diploid nuclear.
##' @return returns a phylo or multiphylo object with branch lengths
##' on the same scale as the species tree.
##' @author Noah Reid, Francois Michonneau
mstree <- function(phy, msdir, nseq, nreps, samplescheme, ploidy = 1) {

    ### gets tip labels for descendents of node in phy
    node.tips <- function(phy, node) {
        n <- length(phy$tip.label)
        if (node <= n) {
            node
        } else {
            d <- phy$edge[which(phy$edge[, 1] == node), 2]
            l <- numeric(length(d))
            for (j in seq_len(d)) {
                if (j <= n) {
                  l[j] <- j
                } else {
                  l[j] <- node.tips(phy, j)
                }
            }
            l
        }
    }

    nspec <- length(phy$tip.label)

    values <- phy$dmv * 2 * ploidy  ###dmv values are now converted to theta=4Nu. dividing branches by theta yields species tree branch lengths in 4N generations, which ms uses for simulation.
    values <- values  ####as on the next line, 'dmv' is Nu haploid alleles, or 2Nu individuals.  Therefore, multiply by 2 to get the standard 4Nu population scaled mutation rate.

    branchscalar <- values[length(values)]
    phy$edge.length <- phy$edge.length/branchscalar  ###transforms branch lengths by theta.  ##BEAST's dmv is equal to Nu where N is the pop size in alleles.  in other words 2Nu in diploid individuals.  so multiply dmv by 2 to get branch lengths in 4N generations for ms.
    values <- values/values[length(values)]  ####dmvs should are now a fraction of the root dmv.
    bt <- sort(branching.times(phy))
    scaleindex <- c(phy$edge[, 2], (nspec + 1))

    #### need to set initial pop sizes -I npops samplescheme[,2] <-n .... -n...>
    #### comlinetree etc...  -n pop scale<size=(scale*No)> -en time pop
    #### scale<newsize=(scale*No)>

    initialpops <- c()
    for (i in 1:nspec) {
        initialpops <- c(initialpops, "-n", i, values[which(phy$edge[, 2] == which(phy$tip.label ==
            samplescheme[i, 1]))])

    }


    comlinetree <- c()
    for (i in 1:length(bt)) {
        child <- sort(phy$edge[phy$edge[, 1] == names(bt[i]), 2])
        tips <- sort(c(sort(node.tips(phy, child[1]))[1], sort(node.tips(phy, child[2]))[1]))
        popi <- which(samplescheme[, 1] == phy$tip.label[tips[2]])
        popj <- which(samplescheme[, 1] == phy$tip.label[tips[1]])
        scale <- values[which(scaleindex == names(bt[i]))]

        comlinetree <- c(comlinetree, "-ej", bt[i], popi, popj, "-en", bt[i], popj,
            scale)
    }
    names(comlinetree) <- NULL

    commandline <- c("cd", msdir, ";", "./ms", nseq, nreps, "-T", "-I", length(phy$tip.label),
        samplescheme[, 2], initialpops, comlinetree, "| grep \\;")
    junk <- system(paste(unlist(commandline), collapse = " "), intern = TRUE)
    trees <- read.tree(text = junk)

    if (nreps > 1) {
        for (i in 1:nreps) {
            trees[[i]]$edge.length <- trees[[i]]$edge.length * branchscalar
        }
    } else {
        trees$edge.length <- trees$edge.length * branchscalar
    }
    return(trees)
    # return(paste(unlist(commandline), collapse = ' '))
}
