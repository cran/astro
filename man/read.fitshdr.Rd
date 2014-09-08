\name{read.fitshdr}
\alias{read.fitshdr}

\title{Read FITS Header}

\description{The mid-level function 'read.fitshdr' allows FITS headers to be read directly into R.}

\usage{read.fitshdr(file, hdu = 1, comments = TRUE, 
             strip = c(" ","'"," "), maxlines = 50000)}

\arguments{
  \item{file}{file name}
  \item{hdu}{header and data unit to be read}
  \item{comments}{output header comments?}
  \item{strip}{lead/trail characters stripped from header 'value' data}
  \item{maxlines}{maximum number of header lines}
}

\details{
The mid-level function 'read.fitshdr' is a wrapper around 'read.fits', and provides a more simplistic output of that routine.
}

\value{
A matrix of data corresponding to the original FITS image header, with two columns ('key' and 'value') and an optional third column ('comment'), should comments = TRUE (default).
}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
