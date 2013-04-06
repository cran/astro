\name{lookback}
\alias{lookback}

\title{Lookback Time}

\description{Calculates the lookback time (light travel time)}

\usage{lookback(z = 1, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Gyr")}

\arguments{
  \item{z}{redshift}
  \item{H}{Hubble constant (km/s/Mpc)}
  \item{M}{Omega M - matter}
  \item{L}{Omega L - energy}
  \item{K}{Omega K - curvature}
  \item{units}{output units [Gyr/s]}
}

\value{
Lookback time in indicated units to the given redshift with the given cosmology.
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

