\name{shadowtext}
\alias{shadowtext}

\title{Add Shadowed Text to a Plot}

\description{Add text with a shadowed background to a plot.}

\usage{shadowtext(x, y = NULL, labels, col = "white", bg = "grey25", 
           theta = (1:8/4)*pi, r1 = 0.06, r2 = 0.04, ...)}

\arguments{
  \item{x,y}{numeric vectors of coordinates where the text labels should be written}
  \item{labels}{a character vector or expression specifying the text to be written}
  \item{col}{foreground text colour}
  \item{bg}{background text colour}
  \item{theta}{position angle offsets for the background shading}
  \item{r1}{shadow width}
  \item{r2}{shadow height}
  \item{...}{additional arguments to be passed to 'text'}
}

\details{
'shadowtext' works by plotting the text initially at positions slightly offset from the central position in a darker background colour, and then over-plots the desired text (twice, to avoid anti-aliasing effects).
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
