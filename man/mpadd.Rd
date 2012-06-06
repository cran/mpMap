\name{mpadd}
\alias{mpadd}
\title{Add markers onto an existing 'mpcross' object}
\usage{
  mpadd(mpcross1, mpcross2, r, theta = 0.15, LOD = 5,
    mapfx = c("haldane", "kosambi"))
}
\arguments{
  \item{mpcross1}{Original object of class \code{mpcross}}

  \item{mpcross2}{Additional object of class
  \code{mpcross}}

  \item{r}{Grid of recombination fractions over which to
  maximize likelihood}

  \item{theta}{Threshold for recombination fractions used
  in constructing linkage groups}

  \item{LOD}{Threshold for LOD scores used in constructing
  linkage groups}

  \item{mapfx}{Map function for converting recombination
  distances to map distances}
}
\value{
  Object of class 'mpcross' which is of the same stage as
  the first object. Recombination fractions and linkage
  groups will be recomputed after adding in the markers in
  the second object if necessary.
}
\description{
  If the 'mpcross' object contains only genotypes and
  pedigree, the markers are added to the genotypes after
  some format checking. If the 'mpcross' object contains
  estimated recombination fractions and LOD scores, the
  markers are added in after estimating recombination
  fractions with the existing markers. If the 'mpcross'
  object contains a map, the markers are added onto the map
  using 3-point mapping.
}
\seealso{
  \code{\link[mpMap]{mpestrf}},
  \code{\link[mpMap]{mpgroup}},
  \code{\link[mpMap]{mpcross}}
}

