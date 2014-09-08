\name{write.fits}
\alias{write.fits}
\alias{.dummy.fits.hdr}
\alias{.check.fits.hdr}
\alias{.make.fits.hdr}

\title{Write FITS Files}

\description{The generic function 'write.fits' allows FITS images (including headers) to be written directly from R.}

\usage{write.fits(x, file = "star.fits", type = "single", hdu = 0)}

\arguments{
  \item{x}{input data (may be a list, see details)}
  \item{file}{file name}
  \item{type}{data format type [single/double/auto]}
  \item{hdu}{write a specific hdu from input list 'x' (0 = all/NULL)}
}

\details{
'write.fits' will write out the data in object x into the named file as a FITS image. x can be in the form of a named list such as that output by 'read.fits', or simply a matrix (or more simply, a vector) of data. 

This function allows custom headers to be used when creating the FITS image. Checks for FITS-critical keywords will be made prior to writing the image, and, if necessary, these keywords will be forcefully changed/removed/added in order to comply with the FITS standard should the header provided violate these rules.

The high-level function 'write.fits' encompasses several low-level functions including '.dummy.fits.hdr', '.check.fits.hdr' and '.make.fits.hdr'. These low-level functions require specific inputs, and must be used in the correct order (particularly in the case of multi-HDU FITS files). For this reason, usage of these low-level functions is not advised in most cases.

Note that writing of FITS binary tables is not yet implemented.
}

\references{
Code and inspiration from the FITSio package <http://cran.r-project.org/web/packages/FITSio/index.html>, written by Andrew Harris <harris at astro.umd.edu>.
}

\author{
Lee Kelvin, Andrew Harris, Aaron Robotham

Maintainer: Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{

require(astro)

# create fake data
dat1 = matrix(rnorm(100*50),100,50)
dat2 = matrix(rnorm(50*25),25,50)

# create multi-HDU FITS image
write.fits(list(dat1,dat2), file="astro.fits")
grep("astro.fits", dir(), value=TRUE)

# read FITS image
x = read.fits("astro.fits")
summary(x)

# show keywords in primary header
x$hdr[[1]][,"key"]

# add keywords into secondary header
write.fitskey(key=c("A","COMMENT"), value=c("B","N/A"), file="astro.fits", 
comment=c("C","astro.fits created by the 'astro' package"), hdu=2)

# print values of 'NAXIS1' and 'A' from secondary header
read.fitskey(c("NAXIS1","A"), "astro.fits", hdu=2)

# create a plot
x = read.fits("astro.fits")
layout(cbind(c(1,2),c(1,3)), widths=c(1,2))
par("mar"=c(3.1,3.1,1.1,1.1))
im1 = x$dat[[1]]
im2 = x$dat[[2]]
image(1:dim(im1)[1], 1:dim(im1)[2], im1, asp=1, xlab="", ylab="")
label("topleft", txt="HDU 1", cex=2, lwd=0)
box()
image(1:dim(im2)[1], 1:dim(im2)[2], im2, asp=1, xlab="", ylab="", 
col=rainbow(1000))
label("topleft", txt="HDU 2", cex=2, lwd=0)
box()
par("mar"=c(3.1,0,1.1,1.1))
aplot(sin, type="n", axes=FALSE, xlab="", ylab="")
hdr = x$hdr[[2]]
ktxt = paste("** astro.fits: HDU 2 Header **\n\nKey\n-----\n",paste(hdr[,"key"],
collapse="\n",sep=""),collapse="",sep="")
vtxt = paste("\n\nValue\n-----\n",paste(hdr[,"value"],collapse="\n",sep=""),
collapse="",sep="")
mtxt = paste("\n\nComment\n-----\n",paste(hdr[,"comment"],collapse="\n",sep=""),
collapse="",sep="")
label("topleft", txt=ktxt, align="left", bty="n")
label("topleft", txt=vtxt, align="left", bty="n", inset=c(1,0.08))
label("topleft", txt=mtxt, align="left", bty="n", inset=c(2,0.08))
label("bottom", txt="note: 'astro.fits' has been automatically deleted", 
bty="n", col="blue", cex=1.5)

unlink("astro.fits")

}

\keyword{data}
