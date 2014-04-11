\name{starbeastPPS-package}
\alias{starbeastPPS-package}
\alias{starbeastPPS}
\docType{package}
\title{
A posterior predictive check for the multispecies coalescent model implemented in *BEAST
}
\description{This package was written to check the fit of the multispecies coalescent model implemented in *BEAST to empirical data. It uses posterior predictive simulation to test whether or not the posterior distribution adequately characterizes the data used to generate it. If it does not, this can be a sign that processes other than stochastic coalescence (such as migration, selection, gene duplication/extinction have shaped individual gene genealogies, or that species assignments have been made incorrectly. If these processes have shaped the data, they can cause *BEAST to return misleading phylogenetic estimates. They are also of fundamental interest to biologists, so poor fit can be a sign that something interesting and worth further investigation is going on in the system. This package was developed along with a meta-analysis of empirical data and was been submitted to Systematic Biology in Dec. 2012 (Reid et al.). For more information on posterior predictive simulation, see Gelman et al.'s book Bayesian Data Analysis. 
}
\details{
\tabular{ll}{
Package: \tab starbeastPPS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-12-02\cr
License: \tab GPLv3\cr
}
}
\author{
Noah M. Reid
Maintainer: Noah M. Reid <noah.reid@gmail.com>
}
\references{
Reid et al. 2012. Poor fit to the multispecies coalescent is widely detectable in empirical data. Submitted to Systematic Biology
Gelman et al. 2003. Bayesian Data Analysis.
}
\examples{
\dontrun{
#data must be read in using read.starbeast()
#this function assumes that the necessary files have the name insertmylocusname.combined.trees, species.combined.trees, or combined.log. do ?read.starbeast for details. 

read.starbeast("path/mydata.XML", "path/directorycontainingtreesandlogfile")->mydata

#a small test dataset, including only 100 samples from the *BEAST posterior (Tamias, from Reid et al. 2012) is included. 

data(tamias)

#the fit of the coalescent genealogies to the species tree can be analyzed by the following. Hudson's ms is required to simulate gene genealogies, and the directory containing the executable must be specified. 

analyze.coalescent(tamias, "path/directorycontainingHudson's_ms")->tamias.coal

#results can be summarized using

summary(tamias.coal)

#and visualized using 

plot(tamias.coal)

#the fit of the sequence data to the gene trees can be analyzed by the following. Rambaut's seq-gen is required to simulate sequence data, and the directory containing the executable must be specified.

analyze.sequences(tamias, "path/directorycontainingSeqGen")->tamias.seq

#results can be summarized using 

summary(tamias.seq)

#and visualized using 

plot(tamias.seq)

}
}
