\name{angdist}
\alias{angdist}

\title{Angular Diameter Distance}

\description{Calculates angular diameter distance}

\usage{angdist(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Mpc")}

\arguments{
  \item{z}{redshift}
  \item{c}{speed of light (m/s)}
  \item{H}{Hubble constant (km/s/Mpc)}
  \item{M}{Omega M - matter}
  \item{L}{Omega L - energy}
  \item{K}{Omega K - curvature}
  \item{units}{output units [Mpc/ly/m]}
}

\value{
Angular diameter distance in indicated units to the given redshift with the given cosmology. Note that this function is only valid for values of K >= 0.
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

