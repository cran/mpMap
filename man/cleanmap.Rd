\name{cleanmap}
\alias{cleanmap}
\title{Clean map for use in QTL mappin}
\usage{
  cleanmap(mpcross, mindist = 1)
}
\arguments{
  \item{mpcross}{Object of class \code{mpcross}}

  \item{mindist}{Minimum distance between markers in cM}
}
\value{
  An mpcross object is returned which has markers removed
  which are within the minimum distance specified. These
  markers do not provide additional information for QTL
  mapping and increase the computational burden.
}
\description{
  Given an mpcross object, this function will remove
  markers which are clustered too tightly together to be
  useful for QTL mapping. Markers will be removed from the
  map, from the data, and from the estimated recombination
  fractions, effectively subsetting the cross down to a
  more spaced out grid of markers.
}
\seealso{
  \code{\link[mpMap]{mpcross}}
}

