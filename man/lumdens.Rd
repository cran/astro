\name{lumdens}
\alias{lumdens}

\title{Luminosity Density}

\description{Calculate the luminosity density from a given set of input (Schechter) luminosity function inputs.}

\usage{lumdens(knee, slope, norm, msun = solar("r"), mag = TRUE, log = FALSE)}

\arguments{
  \item{knee}{The knee(s) of the luminosity distribution}
  \item{slope}{The faint-end slope(s) of the luminosity distribution}
  \item{norm}{The normalisation(s) of the luminosity distribution}
  \item{msun}{The absolute magnitude of the sun}
  \item{mag}{Are the input knee values in magnitudes?}
  \item{log}{Are the input knee values logged?}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
