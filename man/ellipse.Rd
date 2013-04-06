\name{ellipse}
\alias{ellipse}

\title{Calculates Ellipse Coordinates}

\description{Calculates x,y cartesian position coordinates of an ellipse.}

\usage{ellipse(xcen = 0, ycen = 0, a = 10, b = 5, e = 1-b/a, pa = 0)}

\arguments{
  \item{xcen}{origin (x)}
  \item{ycen}{origin (y)}
  \item{a}{semi-major axis}
  \item{b}{semi-minor axis}
  \item{e}{ellipticity (1 - b/a)}
  \item{pa}{position angle (right=0, up=90)}
}

\details{
Note that 'b' is redundant if values of 'e' are given.
}

\value{
  \item{x}{ellipse coordinates along the x-axis}
  \item{y}{ellipse coordinates along the y-axis}
}

\author{
Menaechmus, Euclid, Apollonius, Lee Kelvin

Maintainer: Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
