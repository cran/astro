\name{kronrad}
\alias{kronrad}

\title{Kron Radius}

\description{Calculates the Kron radius based on input Sersic parameters.}

\usage{kronrad(n, r = 1e10, re = 1)}

\arguments{
  \item{n}{Sersic index}
  \item{r}{radius at which to evaluate the Kron radius}
  \item{re}{half-light radius}
}

\details{
The Kron radius is defined as the first moment of the surface brightness light profile. Kron (1980) argued that apertures of multiple values of the Kron radius, Rk (typically, Rk = 2 or Rk = 2.5) would contain a sufficient amount of the total flux of the galaxy to be a useful measure of total flux. This is only true if one is able to evaluate the Kron radius at very large radii away from the galaxy centre (using very deep imaging data). If this is not possible (which typically it is not) the Kron radius will be an underestimate, and consequently, so will the Kron magnitude. This problem becomes most acute for those galaxies with high Sersic indices found within shallow imaging data.
}

\references{Graham A. W., Driver S. P., 2005, PASA, 22, 118}

\author{
Lee Kelvin, <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
