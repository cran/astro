\name{write.fitshdr}
\alias{write.fitshdr}

\title{Write FITS Headers}

\description{The mid-level function 'write.fitshdr' allows FITS headers to be written directly from R.}

\usage{write.fitshdr(hdr, file, hdu = 1)}

\arguments{
  \item{hdr}{input header (see details)}
  \item{file}{file name}
  \item{hdu}{header and data unit to be written to}
}

\details{
Input headers ('hdr') must contain at least two columns, one named 'key' and one named 'value' (a 'comment' column is optional). For an example of the format expected, see 'read.fitshdr' and 'read.fits'.

The mid-level function 'write.fitshdr' is a wrapper around 'read.fits' and 'write.fits', and provides a shortcut to updating the header of an already existent FITS file.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
