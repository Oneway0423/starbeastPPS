\name{mstree}
\alias{mstree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate coalescent genealogies given a species tree using Hudson's ms
}
\description{
takes a species tree of object 'phylo' with a vector giving branch effective population sizes and simulates a coalescent genealogy. this is used internally, but may be useful for some people so it is documented here. 
}
\usage{
mstree(phy, msdir, nseq, nreps, samplescheme, ploidy = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
an ultrametric phylo object with a vector named "dmv" appended giving the effective population sizes. Species tree branch lengths are effectively divided by 2*dmv to put the branch lengths in coalescent units. dmv should therefore be 2Neu (theta/2) if species tree branch lengths are in units of substitutions per site, or 2Ne if branch lengths are in generations. 
}
  \item{msdir}{
the directory containing Hudson's ms. 
}
  \item{nseq}{
the number of alleles to simulate
}
  \item{nreps}{
the number of trees to simulate
}
  \item{samplescheme}{
a data frame with 2 columns and nseq rows. the first column contains the tip labels of the species tree and the second contains the number of alleles to be sampled from the corresponding population. 
}
  \item{ploidy}{
a multiplier for the 'dmv' vector. if ploidy is specified, then dmv=dmv*ploidy. do not specify it unless you wish the dmv vector to be scaled by ploidy. 
}
}
\details{
simulates coalescent genealogies from species trees. effective population sizes on each branch must be specified in a 'dmv' vector appended to the 'phylo' object. 
}
\value{
returns a phylo or multiphylo object with branch lengths on the same scale as the species tree. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Noah Reid
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
