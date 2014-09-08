\name{strip}
\alias{strip}

\title{Strip Leading/Trailing Characters}

\description{Strips leading/trailing characters from the input string. Particularly useful for extracting information from FITS file headers which are embedded in a complex string, for example.}

\usage{strip(x, strip = " ")}

\arguments{
  \item{x}{input string}
  \item{strip}{character to be stripped (may be a vector)}
}

\details{
Characters in the string will be stripped 'outside in', from left-to-right in the order they are given in the argument. See examples below for more detail.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
require("astro")

x = "   'lee'   "

strip(x, strip=" ")
#[1] "'lee'"

strip(x, strip=c(" ","'"))
#[1] "lee"

strip(x, strip=c("'"," "))
#[1] "'lee'"

}

\keyword{data}
