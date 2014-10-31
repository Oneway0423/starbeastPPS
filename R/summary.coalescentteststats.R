summary.coalescentteststats <- function(output) {
    elements <- 2 + length(output$probs$simulated[1, ])
    result <- matrix(nrow = elements, ncol = 2)
    
    probs.g.coef <- (apply(output$probs$simulated, MARGIN = 1, sd)/apply(output$probs$simulated, 
        MARGIN = 1, mean)) - (apply(output$probs$empirical, MARGIN = 1, sd)/apply(output$probs$empirical, 
        MARGIN = 1, mean))
    probs.g <- (output$probs$simulated - output$probs$empirical)
    probs.sum <- rowSums(probs.g)
    
    coal.coef <- (apply(output$coal$simulated, MARGIN = 1, sd)/apply(output$coal$simulated, 
        MARGIN = 1, mean)) - (apply(output$coal$empirical, MARGIN = 1, sd)/apply(output$coal$empirical, 
        MARGIN = 1, mean))
    coal <- (output$coal$simulated - output$coal$empirical)
    coal.sum <- rowSums(coal)
    
    avgs <- mean(probs.sum)
    avgs <- c(avgs, mean(probs.g.coef))
    avgs <- c(avgs, apply(probs.g, MARGIN = 2, mean))
    avgs <- round(avgs, digits = 2)
    
    quants <- quantile(probs.sum, probs = c(0.001, 0.01, 0.025, 0.05, 0.95, 0.975, 
        0.99, 0.999))
    quants <- rbind(quants, quantile(probs.g.coef, probs = c(0.001, 0.01, 0.025, 
        0.05, 0.95, 0.975, 0.99, 0.999)))
    quants <- rbind(quants, t(apply(probs.g, MARGIN = 2, quantile, probs = c(0.001, 
        0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.999))))
    rownames(quants)[1:2] <- c("SUM", "COEF")
    
    # for probabilities, one tailed test.
    values <- quants[, 1:4] > 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 1] <- paste(avgs, values, sep = " ")
    
    avgs <- mean(coal.sum)
    avgs <- c(avgs, mean(coal.coef))
    avgs <- c(avgs, apply(coal, MARGIN = 2, mean))
    avgs <- round(avgs, digits = 2)
    
    quants <- quantile(coal.sum, probs = c(5e-04, 0.005, 0.0125, 0.025, 0.975, 0.9875, 
        0.995, 0.9995))
    quants <- rbind(quants, quantile(coal.coef, probs = c(5e-04, 0.005, 0.0125, 0.025, 
        0.975, 0.9875, 0.995, 0.9995)))
    quants <- rbind(quants, t(apply(coal, MARGIN = 2, quantile, probs = c(5e-04, 
        0.005, 0.0125, 0.025, 0.975, 0.9875, 0.995, 0.9995))))
    rownames(quants)[1:2] <- c("SUM", "COEF")
    
    # for deep coalescences, two tailed test.
    values <- quants[, 1:4] > 0 | quants[, 5:8] < 0
    values <- rowSums(values)
    values[values == 0] <- "n.s."
    values[values == 1] <- "*"  #5%
    values[values == 2] <- "**"  #2.5%
    values[values == 3] <- "***"  #1%
    values[values == 4] <- "****"  #0.1%
    
    result[, 2] <- paste(avgs, values, sep = " ")
    
    rownames(result) <- rownames(quants)
    colnames(result) <- c("gtree.prob", "deep.coal")
    cat("\n\n")
    print(result, quote = FALSE)
    cat("\n\nshown are the means of test statistic distributions and tail area probabilities\n*<0.05, **<0.025, ***<0.01, ****<0.001\ngtree.prob is a one-tailed test. deep.coal is two-tailed.\nall test statistics are calculated as (simulated-empirical) values.\n")
} 
