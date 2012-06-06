\name{computemap}
\alias{computemap}
\title{Computes map distances}
\usage{
  computemap(object, mapfx = c("haldane", "kosambi"),
    missfx = 2)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{mapfx}{Map function to convert recombination
  fractions to cM}

  \item{missfx}{Function to impute missing recombination
  fraction values}
}
\value{
  An mpcross object is returned whose map component has
  been estimated based on the map order and matrix of
  recombination fractions. Missing recombination fractions
  are imputed either by filling in the closest non-missing
  value (missfx=1) or by averaging the distance between
  other nearby markers (missfx=2).
}
\description{
  Given an mpcross object with a map order and matrix of
  recombination fractions, this function will estimate map
  positions
}
\seealso{
  \code{\link[mpMap]{mpcross}}
}

