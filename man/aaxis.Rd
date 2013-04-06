\name{aaxis}
\alias{aaxis}

\title{Add a Scientific Axis to a Plot}

\description{Adds an axis to the current plot.}

\usage{aaxis(side = 1, at = NULL, labels = TRUE, majticks = TRUE, 
      minticks = TRUE, nmaj = NULL, nmin = NULL, unlog = FALSE, 
      format = NULL, digits = NULL, las = NULL, mgp = NULL, 
      tcl = NULL, dexcl = 0.2, ...)}

\arguments{
  \item{side}{an integer specifying which side of the plot the axis is to be drawn on. The axis is placed as follows: 1=below, 2=left, 3=above and 4=right}
  \item{at}{the points at which the major tick-marks are to be drawn}
  \item{labels}{should axis annotation be added}
  \item{majticks}{should major tick marks be included}
  \item{minticks}{should minor tick marks be included}
  \item{nmaj}{number of major tick marks}
  \item{nmin}{number of minor tick marks between major marks}
  \item{unlog}{unlog axis annotation when the data is logged (base 10)}
  \item{format}{format for axes labelling (see 'formatC')}
  \item{digits}{number of digits for axes labels}
  \item{las,mgp,tcl}{standard 'par' plotting parameters}
  \item{dexcl}{distance from major tick marks within which no minor tick labels should be plotted}
  \item{...}{additional arguments to be passed to 'axis'}
}

\details{
The mid-level function 'aaxis' (astro:axis) is a wrapper around the R function 'axis'. It provides significant additional features which trivially allow the creation of figures more suited for a scientific audience. Notably, 'aaxis' allows minor tick-marks and logged axes to be created with the minimum of effort.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
par("mar"=c(5.1,4.1,2.1,4.1))
aplot(1:1000, log10(1:1000), unlog="y", type="l", format="p", side=1:3, col="red", lwd=2)
aaxis(4, nmaj=4, nmin=9)
mtext(bquote(paste(log[10]," y ")), side=4, line=2.5)
label("bottomright", txt="astro:axis (aaxis)", cex=2, lwd=0, bgcol=NULL)
}

\keyword{data}
