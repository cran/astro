\name{scalemark}
\alias{scalemark}

\title{Add a Scalemark to a Plot}

\description{Adds a scalemark to a figure, denoting angular size.}

\usage{scalemark(len = "AUTO", txt = "AUTO", pixsize = NA, 
          col = "white", linecol = col, bg = "grey25", 
          pos = "bottomleft", inset = 0.1, cex = 0.9, lwd = 1.2)}

\arguments{
  \item{len}{angular size to represent ('AUTO' = automatic best fit)}
  \item{txt}{angular size text ('AUTO' = automatic as above)}
  \item{pixsize}{pixel scale (arcsec/pixel)}
  \item{col}{colour}
  \item{linecol}{line colour}
  \item{bg}{shadow colour}
  \item{pos}{position (bottomleft/topleft/topright/bottomright)}
  \item{inset}{inset into plot}
  \item{cex}{text expansion factor}
  \item{lwd}{line widths}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
