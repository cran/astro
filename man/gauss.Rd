\name{gauss}
\alias{gauss}

\title{The Gaussian Function}

\description{This function calculates the Gaussian function for an input data range.}

\usage{gauss(x, mean = 0, sigma = 1, norm = 1, bw = 0.1, ftype = "lin", ...)}

\arguments{
  \item{x}{input values (typically a smooth data range)}
  \item{mean}{the mean(s) of the Gaussian function}
  \item{sigma}{the 1-sigma value}
  \item{norm}{the normalisation(s) of the Gaussian function}
  \item{bw}{integration bin width sizes}
  \item{ftype}{type of input data [lin/log/ln]}
  \item{...}{additional arguments to be passed to 'integrate'}
}

\value{
A vector of length equal to the length of data representing the number density *per dex* at each input data point. Note: to convert the final number densities into their original bin-width (e.g., per 0.5 dex) multiply the output of this function by the bin-width.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
