\name{aqbeta}
\alias{aqbeta}

\title{Astro:Beta Distribution (quantile)}

\description{An astro wrapper around the 'qbeta' function. Useful for calculating robust upper and lower error boundaries for a given data set.}

\usage{aqbeta(k, n, s = c(-1,1), p = NA, corr = TRUE, ...)}

\arguments{
  \item{k}{number of successes}
  \item{n}{total number of trials}
  \item{s}{sigma values required}
  \item{p}{probability values required (see details)}
  \item{corr}{apply a one-sided correction for extreme values ( k=0 / k=n )}
  \item{...}{arguments to be passed to 'qbeta'}
}

\details{
When 'p' is equal to <NA> (default), sigma confidence intervals are calculated according to the value of 's'.  If a value of 'p' is given, probabilities of 'p' are calculated instead.
}

\references{
Cameron, E., 2011, Publications of the Astronomical Society of Australia (PASA),  28, 128
}

\author{Ewan Cameron <dr.ewan.cameron@gmail.com>

Aaron Robotham <aaron.robotham@icrar.org>

Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
