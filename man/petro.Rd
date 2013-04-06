\name{petro}
\alias{petro}

\title{Petrosian Function}

\description{Calculates Petrosian parameters based on input Sersic parameters.}

\usage{petro(mag, n, e = 0, rp = 3, i = 0.5)}

\arguments{
  \item{mag}{total magnitude}
  \item{n}{Sersic index}
  \item{e}{ellipticity (1 - b/a)}
  \item{rp}{integrate out to this many multiples of the Petrosian radius}
  \item{i}{use this value of 1/petrosian index}
}

\details{
The Petrosian function describes the ratio between the average intensity within some projected radius and the intensity at that radius. The value of this ratio is known as the Petrosian index. Typically, a fixed value of 1/Petrosian index, i (usually i = 0.2 or i = 0.5) is chosen in order to define the Petrosian radius. The Petrosian magnitude is then defined as the magnitude lying within a given multiple of the Petrosian radius, Rp (typically, Rp = 2 or Rp = 3 for i = 0.2 and i = 0.5, respectively).
}

\value{
  \item{mag}{magnitude within r}
  \item{magdiff}{difference between total magnitude and magnitude within r}
  \item{mu}{surface brightness at r}
  \item{muavg}{average surface brightness within r}
  \item{inten}{intensity at r}
  \item{lum}{luminosity within r}
  \item{lumtot}{total luminosity}
  \item{lumfrac}{fraction of total luminosity contained within r}
}

\references{Graham A. W., Driver S. P., 2005, PASA, 22, 118}

\author{
Lee Kelvin, <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
