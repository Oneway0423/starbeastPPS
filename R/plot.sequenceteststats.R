plot.sequenceteststats <- function(output, contour = FALSE) {
    devAskNewPage(ask = TRUE)
    
    cat("\n\nAll statistics shown are calculated as (simulated-empirical).\nThe expected mean of each distribution is zero.\n\n")
    elements <- length(output$test.stat$v.sites[1, ])
    samples <- length(output$test.stat$v.sites[, 1])
    
    gcdist <- output$test.stat$gcstat
    sitesdist <- output$test.stat$v.sites
    phy.likdist <- output$test.stat$phy.likelihood
    multi.likdist <- output$test.stat$multi.likelihood
    
    par(mfrow = c(1, 5))
    
    for (i in 1:elements) {
        plot(density(sitesdist[, i]), main = NULL, xlab = "variable sites")
        abline(v = quantile(sitesdist[, i], probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(sitesdist[, i], probs = c(0.005, 0.995)), col = "red", 
            lty = 2)
        abline(v = 0, lwd = 2)
        plot(density(phy.likdist[, i]), main = NULL, xlab = "phylogenetic likelihood")
        abline(v = quantile(phy.likdist[, i], probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(phy.likdist[, i], probs = c(0.005, 0.995)), col = "red", 
            lty = 2)
        abline(v = 0, lwd = 2)
        plot(density(multi.likdist[, i]), main = NULL, xlab = "multinomial likelihood")
        abline(v = quantile(multi.likdist[, i], probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(multi.likdist[, i], probs = c(0.005, 0.995)), col = "red", 
            lty = 2)
        abline(v = 0, lwd = 2)
        plot(density(gcdist[, i]), main = NULL, xlab = "gc statistic")
        abline(v = quantile(gcdist[, i], probs = c(0.025, 0.975)), col = "red")
        abline(v = quantile(gcdist[, i], probs = c(0.005, 0.995)), col = "red", lty = 2)
        abline(v = 0, lwd = 2)
        if (!contour) {
            plot(x = sitesdist[, i], y = gcdist[, i], main = NULL, xlab = "variable sites", 
                ylab = "gc statistic", pch = 21, cex = 0.5, bg = "black")
            abline(v = 0, col = "red")
            abline(h = 0, col = "red")
        } else {
            hdr.boxplot.2d(x = sitesdist[, i], y = gcdist[, i], prob = c(0.001, 0.01, 
                0.05), xlab = "variable sites", ylab = "gc statistic")
            abline(v = 0, col = "red")
            abline(h = 0, col = "red")
            points(x = sitesdist[, i], y = gcdist[, i], cex = 0.5)
        }
    }
    
} 
