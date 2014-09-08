\name{read.fits}
\alias{read.fits}
\alias{.read.fits.hdr}
\alias{.parse.fits.hdr}
\alias{.read.fits.image}
\alias{.read.fits.table}

\title{Read FITS Files}

\description{The generic function 'read.fits' allows FITS binary tables and images (including headers) to be read directly into R.}

\usage{read.fits(file, hdu = 0, comments = TRUE, strip = c(" ","'"," "), 
          maxlines = 50000, xlo = NA, xhi = NA, ylo = NA, yhi = NA)}

\arguments{
  \item{file}{file name}
  \item{hdu}{header and data unit to be read (0 = all)}
  \item{comments}{output header comments?}
  \item{strip}{lead/trail characters stripped from header 'value' data}
  \item{maxlines}{maximum number of header lines}
  \item{xlo}{lower x pixel sub-region (image only)}
  \item{xhi}{upper x pixel sub-region (image only)}
  \item{ylo}{lower y pixel sub-region (image only)}
  \item{yhi}{upper y pixel sub-region (image only)}
}

\details{
The strip argument uses the function 'strip', and removes leading/trailing characters in the order they are given in the argument. The default [c(" ","'"," ")] usually does a good job of removing all trace of FITS header formatting.

The high-level function 'read.fits' encompasses several low-level functions including '.read.fits.hdr', '.parse.fits.hdr', '.read.fits.image' and '.read.fits.table'. These low-level functions require specific inputs, and must be used in the correct order (particularly in the case of multi-HDU FITS files). For this reason, usage of these low-level functions is not advised in most cases.
}

\value{
A list of length two, named '$hdr' and '$dat'. '$hdr' contains the header information for the file, whereas '$dat' contains the imaging or binary table data. Both '$hdr' and '$dat' contain a sub-list for each hdu present in the original file.

Each FITS extension containing a binary table (binary tables are never found in the primary FITS HDU) has additional sub-lists within the '$dat' list, named '$meta' and '$table'. These provide the table meta-data and actual table data itself, respectively.

  \item{$hdr}{Header}
  \item{$dat}{Data Unit}
}

\references{
Code and inspiration from the FITSio package <http://cran.r-project.org/web/packages/FITSio/index.html>, written by Andrew Harris <harris at astro.umd.edu>.
}

\author{
Lee Kelvin, Andrew Harris, Angus Wright

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
ktxt = paste("** astro.fits: HDU 2 Header **\n\nKey\n-----\n",
paste(hdr[,"key"],collapse="\n",sep=""),collapse="",sep="")
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
