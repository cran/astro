\name{put.fitskey}
\alias{put.fitskey}

\title{Put FITS Keyword Value Into Header}

\description{A utility function to allow easy input of a FITS header value into an already loaded FITS header object.}

\usage{put.fitskey(key, value, hdr)}

\arguments{
  \item{key}{header keyword (may be a vector)}
  \item{value}{header value (may be a vector)}
  \item{hdr}{FITS header object}
}

\value{
The updated header object.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
