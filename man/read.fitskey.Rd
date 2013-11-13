\name{read.fitskey}
\alias{read.fitskey}

\title{Read FITS Header Keyword}

\description{The mid-level function 'read.fitskey' allows FITS header keywords to be read directly into R.}

\usage{read.fitskey(key, file, hdu = 1)}

\arguments{
  \item{key}{header keyword (may be a vector)}
  \item{file}{file name}
  \item{hdu}{header and data unit to be read}
}

\details{
The mid-level function 'read.fitskey' is a wrapper around 'read.fits', and provides a more simplistic output of that routine.
}

\value{
A vector of data equal in length to the input key request. NA is returned where no keys have been found.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
