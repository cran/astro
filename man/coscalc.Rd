\name{coscalc}
\alias{coscalc}

\title{Calculates Cosmological Parameters}

\description{Calculates several commonly used cosmological parameters using a given cosmology.}

\usage{coscalc(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, 
        K = 1-M-L, dunit = "Mpc", vunit = "Gpc3", tunit = "Gyr", 
        r = 1, inp = "arcsec", out = "kpc")}

\arguments{
  \item{z}{redshift}
  \item{c}{speed of light (m/s)}
  \item{H}{Hubble constant (km/s/Mpc)}
  \item{M}{Omega M - matter}
  \item{L}{Omega L - energy}
  \item{K}{Omega K - curvature}
  \item{dunit}{outout distance units [Mpc/ly/m]}
  \item{vunit}{output volume units [Gpc3/ly3/m3]}
  \item{tunit}{output time units [Gyr/s]}
  \item{r}{radius/size}
  \item{inp}{input angular size units (see details)}
  \item{out}{output angular size units (see details)}
}

\details{
'coscalc' brings together several commonly used cosmological distance/volume/time measurements as given by several other functions in the 'astro' package. 

Units available for angular size conversion are: 'deg', 'rad' or 'arcsec' <==> 'm', 'pc', 'kpc' or 'Mpc'.
}

\value{
  \item{z}{redshift}
  \item{codist.los}{comoving distance (line-of-sight)}
  \item{codist.trans}{comoving distance (transverse)}
  \item{angdist}{angular diameter distance}
  \item{lumdist}{luminosity distance}
  \item{covol}{comoving volume}
  \item{lookback}{lookback time}
  \item{age}{age of Universe at 'z'}
  \item{angsize}{angular size conversion}
}

\references{Hogg D. W., 1999, ArXiv Astrophysics e-prints, arXiv:astro-ph/9905116}

\author{
Lee Kelvin, Aaron Robotham

Maintainer: Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

