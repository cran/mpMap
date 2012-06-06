\name{mpgroup}
\alias{mpgroup}
\title{Group markers into linkage groups given 2-pt recombination fraction estimates}
\usage{
  mpgroup(mpcross, theta = 0.15, LOD = 5)
}
\arguments{
  \item{mpcross}{an object of class \code{mpcross} which
  includes the component \code{rf} output by
  \code{mp.est.rf}. See \code{\link[mpMap]{mpcross.object}}
  for further details.}

  \item{theta}{Threshold for grouping based on
  recombination fraction value}

  \item{LOD}{Threshold for grouping based on LOD score
  value}
}
\value{
  The original object, with the added component \code{lg}
  which is a list with the following components:
  \item{n.groups}{ The number of linkage groups formed by
  the function} \item{groups}{ Vector labelling each marker
  by assigned linkage group. Missing values mean that the
  marker was linked to more than one group and could not be
  assigned with confidence} \item{LODthresh}{ The LOD
  threshold value used to determine linkage}
  \item{thetathresh}{ The theta threshold value used to
  determine linkage} \item{order}{ A list with a component
  for each constructed linkage group which contains the
  order of the markers within the linkage group}
}
\description{
  Groups markers based on estimated pairwise recombination
  fractions. Linkage groups are built up by adding a marker
  if it satisfies the criteria for linkage with at least
  one other marker in the group.
}
\seealso{
  \code{\link[mpMap]{mpestrf}}

  map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE,
  include.x=FALSE) sim.ped <- sim.mpped(4, 1, 500, 6, 1)
  sim.dat <- sim.mpcross(map=map, pedigree=sim.ped,
  qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6,
  byrow=TRUE), seed=1) dat.rf <- mpestrf(sim.dat) dat.lg <-
  mpgroup(dat.rf)
}

