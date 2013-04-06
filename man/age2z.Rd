\name{age2z}
\alias{age2z}

\title{Redshift of the Universe at a given age}

\description{Calculates the redshift of the Universe at a given age}

\usage{age2z(t = 1, steps = 10, H = 70, M = 0.3, L = 1-M, K = 1-M-L, 
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
This is a very crude inverse variation of the 'age' function. As such, the output result should be used as a guide only.
}

\value{
Redshift of the Universe at a given age with the given cosmology.
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

