\name{gauss.fit}
\alias{gauss.fit}
\alias{.gauss.fit.bin}
\alias{.gauss.fit.dat}
\alias{.gauss.fit.chi}

\title{Fit to the Gaussian Function}

\description{This function calculates a Gaussian function fit to a set of input data.}

\usage{gauss.fit(data, mean = 0, sigma = 1, norm = 1, mean.alt = NA, 
          sigma.alt = NA, norm.alt = NA, meanlo = -Inf, sigmalo = 1e-9, 
          normlo = 0, meanhi = Inf, sigmahi = Inf, normhi = Inf, 
          fixmean = FALSE, fixsigma = FALSE, fixnorm = FALSE, 
          range = NA, lim1 = NA, lim2 = NA, numlim = 1, 
          method = "nlminb", volume = 1, bw = 0.1, ftype = "lin", 
          null = 1E-9, error = "jack", subvol = 10, sampnum = subvol, ...)}

\arguments{
  \item{data}{input data vector}
  \item{mean}{the mean(s) of the Gaussian function}
  \item{sigma}{the 1-sigma value of the Gaussian function}
  \item{norm}{the normalisation(s) of the Gaussian function}
  \item{mean.alt}{alternative mean(s) of the Gaussian function}
  \item{sigma.alt}{alternative sigma(s) of the Gaussian function}
  \item{norm.alt}{alternative normalisation(s) of the Gaussian function}
  \item{meanlo}{a lower bound on the mean parameter}
  \item{sigmalo}{a lower bound on the sigma parameter}
  \item{normlo}{a lower bound on the norm parameter}
  \item{meanhi}{an upper bound on the mean parameter}
  \item{sigmahi}{an upper bound on the sigma parameter}
  \item{normhi}{an upper bound on the norm parameter}
  \item{fixmean}{fix the mean(s)?}
  \item{fixsigma}{fix the first sigma(s)?}
  \item{fixnorm}{fix the first normalisation(s)?}
  \item{range}{data range of interest}
  \item{lim1}{lower data limit for fitting}
  \item{lim2}{upper data limit for fitting}
  \item{numlim}{lower number (per bin) limit for fitting}
  \item{method}{choice of 'nlminb' (recommended) or one of 'optim's minimisation methods (e.g., 'Nelder-Mead')}
  \item{volume}{volume across which the data has been collected}
  \item{bw}{bin width sizes}
  \item{ftype}{type of input data [lin/log/ln]}
  \item{null}{value of null}
  \item{error}{parameter error estimation method [jack/boot]}
  \item{subvol}{number of sub-volumes to split the input data into}
  \item{sampnum}{number of samplings to be made for bootstrapping method}
  \item{...}{additional arguments to be passed to 'integrate'}
}

\value{
  \item{binmid}{bin midpoints}
  \item{num}{number per bin}
  \item{den}{density per bin}
  \item{err}{error per bin}
  \item{errlo}{lower error limit}
  \item{errhi}{upper error limit}
  \item{par}{best fit parameters}
  \item{parlo}{lower error estimates on recovered fit parameters}
  \item{parhi}{upper error estimates on recovered fit parameters}
  \item{chi2}{chi2 value}
  \item{dof}{number of degrees of freedom}
  \item{rchi2}{reduced chi2}
  \item{denlim}{lower density limit used}
  \item{hessian}{optim hessian output}
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

# set up data
data = rnorm(1e4, m=5, s=5)

# fit gaussian
gfit = gauss.fit(data, error="NA")

# print results
cat("mean: ", gfit$par["m1"], "\n", "sigma: ", gfit$par["s1"], 
"norm: ", gfit$par["n1"], "\n", sep="")

# plot data
suppressWarnings(aplot(density(data,weights=rep(1,length(data))), xlim=c(-15,25), 
ylim=c(0,1000), ylab="Number", main=""))
abline(v=5, lty=2)

# plot best-fit results
x = seq(-15,25,len=100)
y = gauss(x, mean=gfit$par["m1"], sigma=gfit$par["s1"], norm=gfit$par["n1"])
lines(x, y, col="red")

# legend
legend("topright", legend=c("data","best fit"), lty=1, col=c("black","red"), 
bty="n")
label("topleft", txt=paste("Input\nmean: 5\nsigma: 5\nnorm: 10000\n\nBest Fit
mean: ", formatC(gfit$par["m1"],format="f",digits=2), "\n", "sigma: ", 
formatC(gfit$par["s1"],format="f",digits=2), "\n", "norm: ", 
formatC(gfit$par["n1"],format="f",digits=2), "\n", sep=""), align="left")

}

\keyword{data}
