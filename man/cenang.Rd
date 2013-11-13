\name{cenang}
\alias{cenang}

\title{Central Angle}

\description{Calculates the angle between two points on a sphere}

\usage{cenang(a1, d1, a2, d2, units = "deg", method = "vincenty")}

\arguments{
  \item{a1}{right ascension of point 1}
  \item{d1}{declination of point 1}
  \item{a2}{right ascension of point 2}
  \item{d2}{declination of point 2}
  \item{units}{input/output units [deg/rad]}
  \item{method}{angle calculation method (see below)}
}

\details{
The central angle describes the angle, from the origin, between two points lying on the surface of a sphere (e.g., the celestial sphere). Three commonly used methods employed in its calculation are: the Spherical Law of Cosines (sloc); the Haversine Formula (haversine), and; the Vincenty Formula (vincenty). The Spherical Law of Cosines suffers from severe rounding errors for small angles (theta < 1E-5). The Haversine Formula generally works well at all angles, but suffers from rounding errors for antipodal points. The Vincenty Formula is accurate for all angles, and is recommended. The three methods may be chosen by specifying `sloc', `haversine' or `vincenty', or by their respective first letters: s, h or v.
}

\value{
The central angle: the angle between two points lying on the surface of a sphere, with reference at the origin/centroid.
}

\references{http://en.wikipedia.org/wiki/Great-circle_distance}

\author{
Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

