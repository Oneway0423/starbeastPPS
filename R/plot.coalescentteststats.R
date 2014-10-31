plot.coalescentteststats <- function(output, show = TRUE) {
    devAskNewPage(ask = TRUE)
    
    elements <- length(output$probs$simulated[1, ])
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
    colnames(coal) <- colnames(probs.g)
    
    data <- list()
    data.probabilities <- list()
    data.deepcoal <- list()
    
    data.probabilities[["sum"]] <- probs.sum
    data.probabilities[["coef"]] <- probs.g.coef
    data.probabilities[["loci"]] <- probs.g
    
    data.deepcoal[["sum"]] <- coal.sum
    data.deepcoal[["coef"]] <- coal.coef
    data.deepcoal[["loci"]] <- coal
    
    data[["probabilities"]] <- data.probabilities
    data[["deepcoal"]] <- data.deepcoal
    
    if (show) {
        cat("\n\nblack line: expected mean given good model fit (zero)\nsolid line: tail area probability = 0.05\ndashed line: tail area probability = 0.01\nif the expected mean falls in the tail, the model is a poor predictor of the data\n\n")
        par(mfrow = c(1, 2))
        ## sum of probabilities
        plot(density(probs.sum), main = "probabilities:sum", xlab = "sum of (simulated-empirical) across loci")
        abline(v = quantile(probs.sum, probs = c(0.05)), col = "red")
        abline(v = quantile(probs.sum, probs = c(0.01)), col = "red", lty = 2)
        abline(v = 0)
        plot(density(coal.sum), main = "deep coal.:sum", xlab = "sum of (simulated-empirical) across loci")
        abline(v = quantile(coal.sum, probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(coal.sum, probs = c(0.005, 0.995)), col = "red", lty = 2)
        abline(v = 0)
        
        ## coefficient of variation
        plot(density(probs.g.coef), main = "probabilities:coef. of var.", xlab = "sum of (simulated-empirical) across loci")
        abline(v = quantile(probs.g.coef, probs = c(0.95)), col = "red")
        abline(v = quantile(probs.g.coef, probs = c(0.99)), col = "red", lty = 2)
        abline(v = 0)
        plot(density(coal.coef), main = "deep coal.:coef. of var.", xlab = "sum of (simulated-empirical) across loci")
        abline(v = quantile(coal.coef, probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(coal.coef, probs = c(0.005, 0.995)), col = "red", lty = 2)
        abline(v = 0)
        
        ## individual loci
        loci <- colnames(probs.g)
        for (i in 1:elements) {
            plot(density(probs.g[, i]), main = paste(loci[i], ": probabilities", 
                sep = ""), xlab = "simulated-empirical values")
            abline(v = quantile(probs.g[, i], probs = c(0.05)), col = "red")
            abline(v = quantile(probs.g[, i], probs = c(0.01)), col = "red", lty = 2)
            abline(v = 0)
            plot(density(coal[, i]), main = paste(loci[i], ": deep coalescences", 
                sep = ""), xlab = "simulated-empirical values")
            abline(v = quantile(coal[, i], probs = c(0.025, 0.975)), col = "red")
            abline(v = quantile(coal[, i], probs = c(0.005, 0.995)), col = "red", 
                lty = 2)
            abline(v = 0)
            
        }
    } else {
        cat("\nreturning test statistic distributions\n\n")
        return(data)
    }
    
} 
