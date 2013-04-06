\name{angsize}
\alias{angsize}

\title{Angular Size Conversion}

\description{Converts between angular and physical size}

\usage{angsize(z = 1, r = 1, inp = "arcsec", out = "kpc", c = 3E8, H = 70, 
        M = 0.3, L = 1-M, K = 1-M-L)}

\arguments{
  \item{z}{redshift}
  \item{r}{radius}
  \item{inp}{input units}
  \item{out}{output units}
  \item{c}{speed of light (m/s)}
  \item{H}{Hubble constant (km/s/Mpc)}
  \item{M}{Omega M - matter}
  \item{L}{Omega L - energy}
  \item{K}{Omega K - curvature}
}

\details{
Units available for conversion are: 'deg', 'rad' or 'arcsec' <==> 'm', 'pc', 'kpc' or 'Mpc'.
}

\value{
Converted size in indicated units to the given redshift with the given cosmology.
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

