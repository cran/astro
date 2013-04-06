\name{lookback2z}
\alias{lookback2z}

\title{Redshift at a given Lookback Time}

\description{Calculates the lookback time (light travel time)}

\usage{lookback2z(t = 1, steps = 10, H = 70, M = 0.3, L = 1-M, K = 1-M-L, 
           units = "Gyr")}

\arguments{
  \item{t}{time}
  \item{steps}{number of accuracy steps}
  \item{H}{Hubble constant (km/s/Mpc)}
  \item{M}{Omega M - matter}
  \item{L}{Omega L - energy}
  \item{K}{Omega K - curvature}
  \item{units}{input units [Gyr/s]}
}

\details{
This is a very crude inverse variation of the 'lookback' function. As such, the output result should be used as a guide only.
}

\value{
Redshift at a given lookback time with the given cosmology.
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

