\name{label}
\alias{label}

\title{Add A Label To A Plot}

\description{Does exactly what it says on the tin.}

\usage{label(pos = "topleft", lab = "label", txt = NULL, inset = 0.1, 
      whitespace = 0.08, col = "black", bgcol = "white", bty = "n", 
      bordercol = "black", lwd = 1, cex = 1, align = "center")}

\arguments{
  \item{pos}{position of label (topleft, top, topright, ...)}
  \item{lab}{contents of the label}
  \item{txt}{see 'lab' (backwards compatibility)}
  \item{inset}{label box inset}
  \item{whitespace}{separation between text and box edge}
  \item{col}{colour}
  \item{bgcol}{background colour}
  \item{bty}{box type (b or n)}
  \item{bordercol}{border colour}
  \item{lwd}{line width of the border}
  \item{cex}{expansion factor}
  \item{align}{text align within box}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
par("mar"=c(5.1,4.1,2.1,2.1))
aplot(sin, xlim=c(0,2*pi), ylim=c(-1.1,1.1), bgcol="lightgoldenrodyellow")
abline(h=0, col="grey75")
label("top", txt="Sine Function", lwd=0, bgcol="grey25", col="white")
label("bottomleft", txt="astro:label (label)", cex=2, lwd=0, bgcol=NULL)
}

\keyword{data}
