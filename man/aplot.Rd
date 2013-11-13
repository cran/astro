\name{aplot}
\alias{aplot}

\title{Scientific X-Y Plotting}

\description{Generic function for plotting of R objects.}

\usage{aplot(x, y = NULL, z = NULL, xlim = NULL, ylim = NULL, zlim = NULL, 
      xlab = NULL, ylab = NULL, zlab = NULL, col = NULL, axes = TRUE, 
      side = 1:4, labels = TRUE, majticks = TRUE, minticks = TRUE, 
      nxmaj = NULL, nymaj = NULL, nxmin = NULL, nymin = NULL, xat = NULL, 
      yat = NULL, log = "", unlog = FALSE, xformat = NULL, yformat = NULL, 
      digits = 0, cex = 1, xlabpos = 1, ylabpos = 2, zcol = NULL, 
      cb = FALSE, cbpos = 4, cbsep = 1.5, cbspan = 2, cbinset = 1, 
      cbx1 = NULL, cbx2 = NULL, cby1 = NULL, cby2 = NULL, cblegend = NULL, 
      cbsteps = 250, las = 0, mgp = c(2.5,0.5,0), tcl = 0.5, dexcl = 0.2, 
      cex.axis = 1, cex.cb = 1, zline = mgp[1]+1, col.axes = "black", 
      col.axis = "black", add = FALSE, type = NULL, bgcol = NULL, ...)}

\arguments{
  \item{x,y,z}{the 'x', 'y' and 'z' arguments provide the x, y and z coordinates for the plot. Supplying the 'z' argument will colour each data point differently according to its z value.}
  \item{xlim,ylim,zlim}{the x, y and z limits of the plot in the form c(lower,upper)}
  \item{xlab,ylab,zlab}{the x, y and z axis labels}
  \item{col}{colour of the data points (will override z)}
  \item{axes}{plot axes}
  \item{side}{sides to plot axes [T/F or 1:4]}
  \item{labels}{sides to plot axes labels [T/F or 1:4]}
  \item{majticks}{plot major tick marks}
  \item{minticks}{plot minor tick marks}
  \item{nxmaj}{number of major tick marks on the x-axis}
  \item{nymaj}{number of major tick marks on the y-axis}
  \item{nxmin}{number of minor tick marks between major ticks (x)}
  \item{nymin}{number of minor tick marks between major ticks (y)}
  \item{xat}{position of x-axes major tick marks}
  \item{yat}{position of y-axes major tick marks}
  \item{log}{logged axes}
  \item{unlog}{unlog axes plotting logged data}
  \item{xformat}{format for x-axes labelling (see 'formatC')}
  \item{yformat}{format for y-axes labelling (see 'formatC')}
  \item{digits}{number of digits for axes labels}
  \item{cex}{expansion factor}
  \item{xlabpos}{label position (x)}
  \item{ylabpos}{label position (y)}
  \item{zcol}{z-axis colour palette}
  \item{cb}{plot colourbar}
  \item{cbpos}{position of colourbar}
  \item{cbsep}{separation of colourbar from plot}
  \item{cbspan}{width/height of colourbar}
  \item{cbinset}{size of inset of colourbar on axis parallel to plotting side}
  \item{cbx1,cbx2,cby1,cby2}{manual colourbar positioning (xlower, xupper, ylower, yupper}
  \item{cblegend}{colourbar legend}
  \item{cbsteps}{colourbar resolution}
  \item{las,mgp,tcl}{standard 'par' plotting parameters}
  \item{dexcl}{distance from major tick marks within which no minor tick labels should be plotted}
  \item{cex.axis}{expansion factor for axis annotation}
  \item{cex.cb}{expansion factor for colourbar annotation}
  \item{zline}{colourbar label line}
  \item{col.axes}{colour of axes}
  \item{col.axis}{colour of axes annotation}
  \item{add}{add the plot to an existing plot}
  \item{type}{plot type}
  \item{bgcol}{background colour}
  \item{...}{additional arguments to be passed to 'plot'}
}

\details{
The top-level function 'aplot' (astro:plot) is a wrapper around the R function 'plot'. It provides significant additional features which trivially allow the creation of figures more suited for a scientific audience. Notably, 'aplot' allows z-axis information to be displayed through the use of colourbars and provides improved axes (including minor-tick marks) through the use of the 'aaxis' function.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\note{
Colourbar features are provided by the 'color.legend' function within the 'plotrix' package (author: Jim Lemon). 'aplot' provides a wrapper around this function, therefore, in order for colourbar features to function correctly, the 'plotrix' package must be installed.
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
# example #1
layout(1)
par("mar"=c(5.1,4.1,2.1,2.1))
par("oma"=c(0,0,0,0))
aplot(sin, xlim=c(0,2*pi), ylim=c(-1.1,1.1), bgcol="lightgoldenrodyellow")
abline(h=0, col="grey75")
label("top", txt="Sine Function", lwd=0, bgcol="grey25", col="white")
label("bottomleft", txt="astro:label (label)", cex=2, lwd=0, bgcol=NULL)

# example #2
layout(1)
par("mar"=c(5.1,4.1,2.1,4.1))
par("oma"=c(0,0,0,0))
aplot(1:1000, log10(1:1000), unlog="y", type="l", yformat="p", side=1:3, 
col="red", lwd=2)
aaxis(4, nmaj=4, nmin=9)
mtext(bquote(paste(log[10]," y ")), side=4, line=2.5)
label("bottomright", txt="astro:axis (aaxis)", cex=2, lwd=0, bgcol=NULL)

# example #3
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

# example #4
layout(1)
par("mar"=c(5.1,4.1,4.1,5.1))
par("oma"=c(0,0,0,0))
aplot(rnorm(1000), rnorm(1000), rnorm(1000), cb=TRUE, zlim=c(-1,1), pch=16, 
xlab="x-axis label", ylab="y-axis label", zlab="z-axis label", bgcol="grey95")
label("topleft", txt="astro:plot (aplot)", cex=2, lwd=0, bgcol=NULL)
}

\keyword{data}
