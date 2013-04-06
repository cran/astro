\name{concen}
\alias{concen}

\title{Central Concentration}

\description{Calculates the central concentration based on input Sersic parameters.}

\usage{concen(a, n)}

\arguments{
  \item{a}{some inner radius (0 < a < 1)}
  \item{n}{Sersic index}
}

\details{
This function calculates the central concentration of a galaxy dependent upon its Sersic index. In order for this calculation to succeed, the outer radial extent of the galaxy must be normalised to 1, and some arbitrary inner radius in the range 0 < a < 1 (for example, a = 1/3) must be chosen.
}

\references{Graham A. W., Driver S. P., 2005, PASA, 22, 118}

\author{
Lee Kelvin, <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
