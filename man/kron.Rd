\name{kron}
\alias{kron}

\title{Kron Function}

\description{Calculates Kron parameters based on input Sersic parameters.}

\usage{kron(mag, n, e = 0, rk = 2.5, r=1e10, re = 1)}

\arguments{
  \item{mag}{total magnitude}
  \item{n}{Sersic index}
  \item{e}{ellipticity (1 - b/a)}
  \item{rk}{integrate out to this many multiples of the Kron radius}
  \item{r}{radius at which to evaluate the Kron radius}
  \item{re}{half-light radius}
}

\details{
The Kron radius is defined as the first moment of the surface brightness light profile. Kron (1980) argued that apertures of multiple values of the Kron radius, Rk (typically, Rk = 2 or Rk = 2.5) would contain a sufficient amount of the total flux of the galaxy to be a useful measure of total flux. This is only true if one is able to evaluate the Kron radius at very large radii away from the galaxy centre (using very deep imaging data). If this is not possible (which typically it is not) the Kron radius will be an underestimate, and consequently, so will the Kron magnitude. This problem becomes most acute for those galaxies with high Sersic indices found within shallow imaging data.
}

\value{
  \item{mag}{magnitude within r}
  \item{magdiff}{difference between total magnitude and magnitude within r}
  \item{mu}{surface brightness at r}
  \item{muavg}{average surface brightness within r}
  \item{inten}{intensity at r}
  \item{lum}{luminosity within r}
  \item{lumtot}{total luminosity}
  \item{lumfrac}{fraction of total luminosity contained within r}
}

\references{Graham A. W., Driver S. P., 2005, PASA, 22, 118}

\author{
Lee Kelvin, <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
