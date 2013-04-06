\name{petroindex}
\alias{petroindex}

\title{Petrosian Index}

\description{Calculates the Petrosian index based on input Sersic parameters.}

\usage{petroindex(r, n, re = 1)}

\arguments{
  \item{r}{radius}
  \item{n}{Sersic index}
  \item{re}{half-light radius}
}

\details{
The Petrosian function describes the ratio between the average intensity within some projected radius and the intensity at that radius. The value of this ratio is known as the Petrosian index. Typically, a fixed value of 1/Petrosian index, i (usually i = 0.2 or i = 0.5) is chosen in order to define the Petrosian radius. The Petrosian magnitude is then defined as the magnitude lying within a given multiple of the Petrosian radius, Rp (typically, Rp = 2 or Rp = 3 for i = 0.2 and i = 0.5, respectively).
}

\references{Graham A. W., Driver S. P., 2005, PASA, 22, 118}

\author{
Lee Kelvin, <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
