\name{astro}
\alias{astro}
\alias{astro-package}
\docType{package}

\title{Astronomy Functions, Tools and Routines}

\description{A collection of commonly used astronomy functions, tools and routines.}

\details{
The astro package provides a series of functions, tools and routines in everyday use within astronomy.

Broadly speaking, one may group these functions into 4 main areas, namely; cosmology, FITS file manipulation, the Sersic function and general (plotting and scripting) tools. A summary of these sub-packages and their functions is as follows:

\subsection{Cosmology}{

Cosmology Calculator: \code{\link{coscalc}},

Comoving distance (line of sight): \code{\link{comovdist.los}},

Comoving distance (transverse): \code{\link{comovdist.trans}},

Angular diameter distance: \code{\link{angdist}},

Luminosity distance: \code{\link{lumdist}},

Comoving volume: \code{\link{comovvol}},

Lookback time: \code{\link{lookback}},

Redshift at a given Lookback time: \code{\link{lookback2z}},

Age of the Universe: \code{\link{age}},

Redshift at a given age: \code{\link{age2z}},

Angular size conversion: \code{\link{angsize}}
    
}

\subsection{FITS}{

Read FITS files (images/tables & headers): \code{\link{read.fits}},

Read FITS header: \code{\link{read.fitshdr}},

Read FITS key: \code{\link{read.fitskey}},

Read FITS image: \code{\link{read.fitsim}},

Read FITS table: \code{\link{read.fitstab}},

Write FITS files: \code{\link{write.fits}},

Write FITS header: \code{\link{write.fitshdr}},

Write FITS key: \code{\link{write.fitskey}},

Plot FITS images: \code{\link{plotfits}},

Strip leading/trailing characters: \code{\link{strip}}

}

\subsection{Sersic}{

The Sersic function: \code{\link{sersic}},

Convert between Sersic radii: \code{\link{convrad}},

Convert between Sersic surface brightnesses: \code{\link{convmu}},

Convert half-light radius to scalelength: \code{\link{re2h}},

Convert scalelength to half-light radius: \code{\link{h2re}},

The incomplete gamma function: \code{\link{igamma}},

Central Concentration: \code{\link{concen}},

The Petrosian function: \code{\link{petro}},

Petrosian index: \code{\link{petroindex}},

Petrosian radius: \code{\link{petrorad}},

The Kron function: \code{\link{kron}},

Kron radius: \code{\link{kronrad}}

}

\subsection{General}{

Astro plot: \code{\link{aplot}},

Astro axis: \code{\link{aaxis}},

Astro colourbar: \code{\link{acb}},

Astro box: \code{\link{abox}},

Astro kernel density estimator: \code{\link{akde2d}},

Astro qbeta function: \code{\link{aqbeta}},

Astro alpha-transparency colour conversion: \code{\link{acol}},

The Gaussian Function: \code{\link{gauss}},

Fit to the Gaussian Function: \code{\link{gauss.fit}},

The Schechter Function: \code{\link{schechter}},

Fit to the Schechter Function: \code{\link{schechter.fit}},

Error ellipses for a Schechter fit: \code{\link{schechter.ellipse}},

Absolute Magnitude of the Sun: \code{\link{solar}},

The Luminosity Density: \code{\link{lumdens}},

Central Angle: \code{\link{cenang}},

Ellipse: \code{\link{ellipse}},

Plot labelling: \code{\link{label}},

Plot shading: \code{\link{shade}},

Human-readable time strings: \code{\link{nicetime}},

For-loop ETA: \code{\link{eta}},

Sigma-clipped mean: \code{\link{scmean}},

Shadowed text: \code{\link{shadowtext}},

Cardinal points: \code{\link{cardinal}},

Scale marks: \code{\link{scalemark}}

}

}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\keyword{package}
\keyword{data}

