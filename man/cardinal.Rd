\name{cardinal}
\alias{cardinal}

\title{Add Cardinal Points to a Plot}

\description{Adds the four cardinal points to a figure, marking north and east.}

\usage{cardinal(rotate = 0, pos = "bottomright", inset = 0.5, 
         len = 2, gap = 1, col = "white", linecol = col, 
         bg = "grey25", lwd = 1, invert = TRUE, south = TRUE, 
         west = TRUE, textrot = TRUE, ...)}

\arguments{
  \item{rotate}{rotate cardinal cross clockwise from 'up' (degrees)}
  \item{pos}{position of cross (bottomleft/topleft/topright/bottomright)}
  \item{inset}{inset cross from edge}
  \item{len}{length of cross arrows}
  \item{gap}{size of gap between arrowhead and text}
  \item{col}{colour of arrows/text}
  \item{linecol}{colour of cardinal cross}
  \item{bg}{shading colour}
  \item{lwd}{line width}
  \item{invert}{invert east-west}
  \item{south}{show the southern prong}
  \item{west}{show the western prong}
  \item{textrot}{rotate the text along with the cross}
  \item{...}{arguments to be passed to 'shadowtext'}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
