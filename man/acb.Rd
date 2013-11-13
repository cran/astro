\name{acb}
\alias{acb}

\title{Add a Colourbar to a Plot}

\description{Adds a colourbar (colour legend) to a plot.}

\usage{acb(zlim = NULL, zlab = NULL, unlog = FALSE, cex = 1, 
    zcol = NULL, cbpos = 4, cbsep = 0, cbspan = 0, cbinset = 1, 
    cbx1 = NULL, cbx2 = NULL, cby1 = NULL, cby2 = NULL, 
    cblegend = NULL, cex.cb = 1, zline = 2.5, ...)}

\arguments{
  \item{zlim}{lower/upper limits of the z-axis}
  \item{zlab}{colourbar label}
  \item{unlog}{unlog logged data}
  \item{cex}{expansion factor}
  \item{zcol}{colourbar colour palette}
  \item{cbpos}{colourbar position (1/2/3/4)}
  \item{cbsep}{separation of colourbar from main plot}
  \item{cbspan}{width/height of colourbar}
  \item{cbinset}{size of colourbar inset parallel to plotting axis}
  \item{cbx1,cbx2,cby1,cby2}{manual placement of colourbar (xlower, xupper, ylower, yupper)}
  \item{cblegend}{colourbar annotation}
  \item{cex.cb}{colourbar expansion factor}
  \item{zline}{colourbar label line}
  \item{...}{additional arguments to be passed to 'color.legend'}
}

\details{
The mid-level function 'acb' (astro:colourbar) is a wrapper around the 'plotrix' function 'color.legend'. It allows trivial creation and placement of colourbars, able to plug-in more readily with figures created with 'aplot'. The main advantage for using 'acb' over that within 'aplot' is for those cases when multi-panelled figures are being created, and the colourbar serves several sub-plots simultaneously.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
layout(cbind(c(1,2),c(3,3)), widths=c(5,1))
par("mar"=c(0,0,0,1))
par("oma"=c(3.1,3.1,3.1,2.1))
aplot(rnorm(1000), rnorm(1000), rnorm(1000), pch=17, zlim=c(-1,1), xlim=c(-3,3), 
ylim=c(-3,3), labels=2:3)
grid()
label("topleft", txt="astro:colourbar (acb)", cex=2, lwd=0, bgcol=NULL)
aplot(rnorm(1000), rnorm(1000), rnorm(1000), pch=16, zlim=c(-1,1), xlim=c(-3,3), 
ylim=c(-3,3), labels=1:2)
grid()
acb(zlim=c(-1,1), zlab="z-axis label")
}

\keyword{data}
