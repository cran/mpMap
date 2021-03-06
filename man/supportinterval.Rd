\name{supportinterval}
\alias{supportinterval}
\title{Calculate support interval for detected QTL}
\usage{
  supportinterval(x, chr, lodsupport = 1)
}
\arguments{
  \item{x}{Object of class \code{mpqtl}}

  \item{chr}{Selected chromosomes}

  \item{lodsupport}{Size of support interval; default is 1
  LOD}
}
\value{
  A list with two components: a matrix containing lower and
  upper bounds for the support intervals for each QTL, and
  the positions of QTL on each chromosome.
}
\description{
  Calculates support interval for QTL based on Wald profile
  and QTL position
}
\details{
  Computes the x-LOD support interval as the region
  surrounding a QTL peak in which the Wald profile exceeds
  the equivalent of x LOD less than the peak value.
}
\examples{
map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
si <- supportinterval(mpq.dat)
}
\seealso{
  \code{\link[mpMap]{plot.mpqtl}}
}

