\name{subset.mpprob}
\alias{subset.mpprob}
\title{Subset mpprob object}
\usage{
  \method{subset}{mpprob} (x, chr = NULL, markers = NULL,
    lines = NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{mpprob}}

  \item{chr}{Selected chromosomes TO KEEP. Default is all}

  \item{markers}{Selected markers TO KEEP. Default is all}

  \item{lines}{Selected lines TO KEEP. Default is all}

  \item{...}{Additional arguments}
}
\value{
  The original object with chromosomes/lines removed which
  are not listed in the arguments.
}
\description{
  Reduces an mpprob object down to a specified set of
  chromosomes, markers and/or lines
}
\note{
  Chromosomes can be input either as the character names of
  chromosomes or the index of the chromosomes in the map.
  Markers can be input as either character values matching
  the colnames of x$finals, or indices of columns in that
  matrix. Note that if markers are removed, the founder
  probabilities will be recomputed for the new map with
  previous settings for mpprob. Lines can be input as
  either character values (matching the rownames of
  x$finals) or indices of rows in that matrix.
}
\seealso{
  \code{\link[mpMap]{mpprob}}
}

