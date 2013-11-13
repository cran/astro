\name{schechter}
\alias{schechter}

\title{The Schechter Function}

\description{This function calculates the single or double Schechter function, most usually associated with the luminosity function.}

\usage{schechter(x, knee, slope, norm, bw = 0.1, mag = FALSE, log = FALSE, ...)}

\arguments{
  \item{x}{input values}
  \item{knee}{the knees(s) of the Schechter function (L_star/M_star)}
  \item{slope}{the slope(s) of the Schechter function (alpha)}
  \item{norm}{the normalisation(s) of the Schechter function (phi_star)}
  \item{bw}{integration bin width sizes}
  \item{mag}{are the input data magnitudes?}
  \item{log}{are the input data logged?}
  \item{...}{additional arguments to be passed to 'integrate'}
}

\value{
A vector of length equal to the length of data representing the number density *per dex* at each input data point. Note: to convert the final number densities into their original bin-width (e.g., per 0.5 dex) multiply the output of this function by the bin-width.
}

\references{
Blanton M. R., Lupton R. H., Schlegel D. J., Strauss M. A., Brinkmann J., Fukugita M., Loveday J., 2005, ApJ, 631, 208

Driver S. P., Popescu C. C., Tuffs R. J., Graham A. W., Liske J., Baldry I., 2008, ApJ, 678, L101

Baldry I. K. Driver S. P., Loveday J., et al., 2012, MNRAS, 421, 621
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{

#
#
# Driver et al. 2008, Figure 2 (magnitude)
#
#

# setup input magnitudes
mag = seq(-24,-17,len=100)

# calculate number densities
num = schechter(mag, bw = 0.5, mag = TRUE, log = FALSE, knee = -21.32, 
slope = -1.32, norm = 4.8E-3)

# plot
aplot(mag, log10(0.5*num), las=1, ylim=c(-6,-2), type="l", xaxs="i", yaxs="i", 
nxmin=1, xlab="Magnitude / mag", 
ylab=bquote(paste("log ", phi, " [ 0.5mag ]", " "^-1)), 
main="Driver et al. 2008, Figure 2")
label("bottomright", lab=bquote(paste("M* = -21.32     ", 
alpha, " = -1.32     ", phi, "* = 4.8x", 10^{-3})), inset=0.1)

#
#
# Baldry et al. 2012, Figure 13 (stellar mass)
#
#

# setup input masses
mass = seq(7,11.6,len=100)

# calculate number densities
num = schechter(mass, bw = 0.1, mag = FALSE, log = TRUE, knee = 10.66, 
slope = c(-0.35,-1.47), norm = c(3.96E-3,0.79E-3))

# plot
aplot(log10(mass), num, las=1, ylim=c(1e-5,2e-1), type="l", log="y", xaxs="i", 
yaxs="i", xlab="log (M / Msun)", 
ylab=bquote(paste("number density (", dex^{-1}, " ", , Mpc^{-3}, ")")), 
yformat="p", main="Baldry et al. 2012, Figure 13")
label("bottomleft", 
lab=bquote(paste("log M* = 10.66   ", alpha[1], " = -0.35   ", 
phi, "*", ""[1], " = 3.96x", 10^{-3}, "   ", alpha[2], " = -1.47   ", 
phi, "*", ""[2], " = 0.79x", 10^{-3})), inset=0.1)

}

\keyword{data}
