\name{plotlink.map}
\alias{plotlink.map}
\title{Plots linkage maps}
\usage{
  plotlink.map(object, chr, max.dist, marker.names = TRUE,
    tick = FALSE, squash = TRUE, colqtl = "red", ...)
}
\arguments{
  \item{object}{Either \code{mpcross} or \code{map} object}

  \item{chr}{Chromosomes to plot}

  \item{max.dist}{Plotting paramters. See
  \code{\link[wgaim]{link.map.cross}}}

  \item{marker.names}{Whether to plot marker names}

  \item{tick}{Plotting parameters. See
  \code{\link[wgaim]{link.map.cross}}}

  \item{squash}{Plotting parameters. See
  \code{\link[wgaim]{link.map.cross}}}

  \item{colqtl}{Color to plot QTL regions. See
  \code{\link[mpMap]{qtlmap}}}

  \item{...}{Additional arguments}
}
\value{
  Modification of link.map.cross function from wgaim to
  allow more general input objects and to highlight regions
  around QTL. If any markers are labelled "QTLx" then they
  will be plotted in a different color.
}
\description{
  Plot linkage map (either as input object or as stored in
  mpcross object). Can also highlight QTL regions when used
  with qtlmap function.
}
\seealso{
  \code{\link[wgaim]{link.map.cross}},
  \code{\link[mpMap]{qtlmap}}
}

