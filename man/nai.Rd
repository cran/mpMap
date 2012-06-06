\name{nai}
\alias{nai}
\title{Count how many generations of advanced intercross are in a pedigree}
\usage{
  nai(pedigree)
}
\arguments{
  \item{pedigree}{Pedigree for a multi-parent cross. Can be
  generated using \code{\link[mpMap]{sim.mpped}}}
}
\value{
  Integer - number of generations of advanced intercrossing
  after mixing stage but before selfing.
}
\description{
  Counts the number of generations of breeding preceding
  selfing and subtracts off the number necessary to
  minimally mix the founders' genomes
}
\examples{
sim.map <- list(Chr1=seq(0,100,10))
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
nai(sim.ped)
sim.ped <- sim.mpped(4, 1, 500, 6, 1, 5)
nai(sim.ped)
}
\seealso{
  \code{\link[mpMap]{sim.mpped}}
}

