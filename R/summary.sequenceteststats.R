summary.sequenceteststats <- function(output) {
    
    cat("\n\nAll statistics shown are calculated as (simulated-empirical).\nThe expected value is zero.\nShown are the means of test statistic distributions and tail area probabilities\n\n")
    elements <- length(output$test.stat$v.sites[1, ])
    samples <- length(output$test.stat$v.sites[, 1])
    
    result <- matrix(nrow = elements, ncol = 4)
    
    
    gcdist <- output$test.stat$gcstat
    sitesdist <- output$test.stat$v.sites
    phy.likdist <- output$test.stat$phy.likelihood
    multi.likdist <- output$test.stat$multi.likelihood
    
    # columns: sites, phy.lik, multi.lik, gc stat treating all as two-tailed tests.
    
    avg.sites <- apply(sitesdist, MARGIN = 2, mean)
    avg.phy <- apply(phy.likdist, MARGIN = 2, mean)
    avg.multi <- apply(multi.likdist, MARGIN = 2, mean)
    avg.gc <- apply(gcdist, MARGIN = 2, mean)
    
    quants.sites <- t(apply(sitesdist, MARGIN = 2, quantile, probs = c(5e-04, 0.005, 
        0.0125, 0.025, 0.975, 0.9875, 0.995, 0.9995)))
    quants.phy <- t(apply(phy.likdist, MARGIN = 2, quantile, probs = c(5e-04, 0.005, 
        0.0125, 0.025, 0.975, 0.9875, 0.995, 0.9995)))
    quants.multi <- t(apply(multi.likdist, MARGIN = 2, quantile, probs = c(5e-04, 
        0.005, 0.0125, 0.025, 0.975, 0.9875, 0.995, 0.9995)))
    quants.gc <- t(apply(gcdist, MARGIN = 2, quantile, probs = c(5e-04, 0.005, 0.0125, 
        0.025, 0.975, 0.9875, 0.995, 0.9995)))
    
    values <- quants.sites[, 1:4] > 0 | quants.sites[, 5:8] < 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 1] <- paste(round(avg.sites, digits = 2), values, sep = " ")
    
    values <- quants.phy[, 1:4] > 0 | quants.phy[, 5:8] < 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 2] <- paste(round(avg.phy, digits = 2), values, sep = " ")
    
    values <- quants.multi[, 1:4] > 0 | quants.multi[, 5:8] < 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 3] <- paste(round(avg.multi, digits = 2), values, sep = " ")
    
    values <- quants.gc[, 1:4] > 0 | quants.gc[, 5:8] < 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 4] <- paste(round(avg.gc, digits = 2), values, sep = " ")
    
    # rownames (loci) not yet part of seq output
    colnames(result) <- c("var. sites", "phy.likelihood", "multin.likelihood", "GC stat")
    cat("\n")
    print(result, quote = FALSE)
    cat("\n\n*<0.05, **<0.025, ***<0.01, ****<0.001\n")
    
    
} 
