\name{read.fitstab}
\alias{read.fitstab}

\title{Read FITS Table}

\description{The mid-level function 'read.fitstab' allows FITS binary tables to be read directly into R.}

\usage{read.fitstab(file, hdu = 2, strip = c(" ","'"," "), maxlines = 50000)}

\arguments{
  \item{file}{file name}
  \item{hdu}{header and data unit to be read}
  \item{strip}{lead/trail characters stripped from header 'value' data}
  \item{maxlines}{maximum number of header lines}
}

\details{
The mid-level function 'read.fitstab' is a wrapper around 'read.fits', and provides a more simplistic output of that routine. Note that the FITS table is usually stored in the 2nd HDU of a FITS file, hence the default setting of hdu for this function.
}

\value{
A matrix of data corresponding to the original FITS table, with column names labelled appropriately.
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
