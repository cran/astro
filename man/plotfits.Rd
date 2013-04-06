\name{plotfits}
\alias{plotfits}

\title{Plot FITS images}

\description{This function allows a FITS image to be converted into a PNG/X11 image, or output as a data array for use within a function. 3 colour RGB images may also be constructed by using multiple file inputs.}

\usage{plotfits(input, hdu = 1, func = "atan", slide = c(0,0,0), 
         scale = c(500,300,100), locut = 0, hicut = pi/2, invert = FALSE, 
         method = 1, type = "x11", width = 5, height = 5, units = "in", 
         res = 300, cen = c(NA,NA), xdim = NA, ydim = NA, 
         file = "image.png")}

\arguments{
  \item{input}{input file(s)}
  \item{hdu}{input hdu}
  \item{func}{scaling function for plot [lin/log/atan]}
  \item{slide}{offset counts from the origin by a given value}
  \item{scale}{scaling value for each image}
  \item{locut}{lower cut (black)}
  \item{hicut}{upper cut (white)}
  \item{invert}{invert the output greyscale/colours}
  \item{method}{counts -> image conversion method [1/2]}
  \item{type}{output type [dat/png/eps/pdf/x11 or bitmap types]}
  \item{width}{output width}
  \item{height}{output height}
  \item{units}{output units}
  \item{res}{output resolution}
  \item{cen}{centre of image}
  \item{xdim}{x width (pixels)}
  \item{ydim}{y width (pixels)}
  \item{file}{output file name}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
