\name{get.fitskey}
\alias{get.fitskey}

\title{Get FITS Keyword Value From Header}

\description{A utility function to allow easy extraction of a FITS header value from an already loaded FITS header object.}

\usage{get.fitskey(key, hdr)}

\arguments{
  \item{key}{header keyword (may be a vector)}
  \item{hdr}{FITS header object}
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
