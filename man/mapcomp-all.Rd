\name{mapcomp}
\alias{mapcomp}
\alias{plot.mapcomp}
\alias{summary.mapcomp}
\title{Functions for comparison of two map orders}
\usage{
mapcomp(object1, object2) 
\method{plot}{mapcomp}(x, ...)
\method{summary}{mapcomp}(object, ...)
}
\arguments{
  \item{object1}{Object inheriting class \code{mpcross} or
  class \code{map}}

  \item{object2}{Object inheriting class \code{mpcross} or
  class \code{map}}

  \item{object}{Object of class \code{mapcomp} for
  summarizing}

  \item{x}{Object of class \code{mapcomp}}

  \item{...}{Additional arguments}
}
\value{
  An object of class \code{mapcomp} with components:
  \item{commonmrk}{A matrix containing 5 columns with the
  names of all common markers for the two maps (mname), the
  chromosome mapped to in the first map (chr1), the
  position mapped to in the first map (pos1), the
  chromosome mapped to in the second map (chr2) and the
  position mapped to in the second map (pos2)}
  \item{samechr}{A matrix containing 5 columns as above.
  Differs from commonmrk in that duplicated markers in
  either map will have been removed, so all markers map to
  exactly one chromosome} \item{map1}{The first map -
  either the originally input object1, or object1$map if it
  inherits class \code{mpcross}} \item{map2}{The first map
  - either the originally input object2, or object2$map if
  it inherits class \code{mpcross}} \item{correlations}{The
  correlation between positions in map1 and map2 for each
  chromosome} \item{dup1}{The names of markers duplicated
  in map1} \item{dup2}{The names of markers duplicated in
  map2}

  Plot produces for a comparison for each chromosome of
  positions of markers which are mapped to that chromosome
  in both maps

  Summary function returns printed summary including -
  number of markers in each map; number of markers common
  to both maps; number of duplicated markers in each map;
  number of markers mapped to different chromosomes;
  correlations between positions on each chromosome.
}
\description{
  Takes in two maps with the aim of comparing the position
  of common markers. Creates a mapcomp object for plotting
  and summary.
}
\examples{
map1 <- sim.map(len=rep(100, 4), n.mar=51, include.x=FALSE)
map2 <- sim.map(len=rep(100, 4), n.mar=51, include.x=FALSE)
mc <- mapcomp(map1, map2)
summary(mc)
plot(mc)
}

