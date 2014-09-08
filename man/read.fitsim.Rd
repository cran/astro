\name{read.fitsim}
\alias{read.fitsim}

\title{Read FITS Image}

\description{The mid-level function 'read.fitsim' allows FITS images to be read directly into R.}

\usage{read.fitsim(file, hdu = 1, maxlines = 50000, xlo = NA, xhi = NA, 
            ylo = NA, yhi = NA)}

\arguments{
  \item{file}{file name}
  \item{hdu}{header and data unit to be read}
  \item{maxlines}{maximum number of header lines}
  \item{xlo}{lower x pixel sub-region (image only)}
  \item{xhi}{upper x pixel sub-region (image only)}
  \item{ylo}{lower y pixel sub-region (image only)}
  \item{yhi}{upper y pixel sub-region (image only)}
}

\details{
The mid-level function 'read.fitsim' is a wrapper around 'read.fits', and provides a more simplistic output of that routine.
}

\value{
A matrix of data corresponding to the original FITS image.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
