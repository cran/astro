\name{write.fitskey}
\alias{write.fitskey}

\title{Write FITS Header Keywords}

\description{The mid-level function 'write.fitskey' allows FITS header keywords to be written directly from R.}

\usage{write.fitskey(key, value, file, comment = "", hdu = 1)}

\arguments{
  \item{key}{header keyword (may be a vector)}
  \item{value}{header value (may be a vector)}
  \item{file}{file name}
  \item{comment}{header comments (may be a vector)}
  \item{hdu}{header and data unit to be written to}
}

\details{
The mid-level function 'write.fitskey' is a wrapper around 'read.fits' and 'write.fits', and provides a shortcut to updating individual keywords of an already existent FITS file.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
