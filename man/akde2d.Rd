\name{akde2d}
\alias{akde2d}

\title{Two-Dimensional Kernel Density Estimation (Astro)}

\description{Two-dimensional kernel density estimation with an axis-aligned bivariate normal kernel, evaluated on a square grid. The function 'akde2d' (astro:kde2d) is a wrapper around the 'kde2d' function within the 'MASS' package. This function adds additional output to that function.}

\usage{akde2d(x, y, n = 25, lims = c(range(x), range(y)), levels = seq(0,0.9,by=0.1), 
       ...)}

\arguments{
  \item{x}{x coordinate of data}
  \item{y}{y coordinate of data}
  \item{n}{number of grid points in each direction. Can be scalar or a length-2 integer vector}
  \item{lims}{the limits of the rectangle covered by the grid as 'c(xl, xu, yl, yu)'}
  \item{levels}{output levels containing the given percentiles of the data}
  \item{...}{arguments to be passed to 'kde2d' in the 'MASS' package}
}

\details{
Characters in the string will be stripped 'outside in', from left-to-right in the order they are given in the argument. See examples below for more detail.
}

\value{
A list of four components.
  
  \item{x,y}{the x and y coordinates of the grid points, vectors of length 'n'}
  \item{z}{an 'n[1]' by 'n[2]' matrix of the estimated density: rows correspond to the value of 'x', columns to the value of 'y'}
  \item{l}{percentile levels containing given fractions of the data}
  
}

\references{
Venables, W. N. and Ripley, B. D. (2002) _Modern Applied Statistics with S._ Fourth edition.  Springer.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
# See 'kde2d' for further examples.
}

\keyword{data}
