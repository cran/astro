\name{schechter.fit}
\alias{schechter.fit}
\alias{.schechter.fit.dat}
\alias{.schechter.fit.chi}

\title{Fit to the Schechter Function}

\description{This function calculates a Schechter function fit to a set of input data.}

\usage{schechter.fit(data, vmax = NA, knee, slope, norm, knee.alt = NA, 
              slope.alt = NA, norm.alt = NA, kneelo = -Inf, 
              slopelo = -Inf, normlo = 0, kneehi = Inf, slopehi = Inf, 
              normhi = Inf, fixk1 = FALSE, fixs1 = FALSE, fixn1 = FALSE, 
              fixk2 = FALSE, fixs2 = FALSE, fixn2 = FALSE, 
              range = range(data), lim1 = NA, lim2 = NA, numlim = 1, 
              method = "nlminb", volume = max(vmax), bw = 0.1, 
              mag = FALSE, log = FALSE, null = 1E-9, error = "jack", 
              subvol = 10, sampnum = subvol, msun = solar("r"))}

\arguments{
  \item{data}{input data vector}
  \item{vmax}{vector of maximum comoving volumes within which object could lie}
  \item{knee}{the knees(s) of the Schechter function (L_star/M_star)}
  \item{slope}{the slope(s) of the Schechter function (alpha)}
  \item{norm}{the normalisation(s) of the Schechter function (phi_star)}
  \item{knee.alt}{alternative knees(s) of the Schechter function (L_star/M_star)}
  \item{slope.alt}{alternative slope(s) of the Schechter function (alpha)}
  \item{norm.alt}{alternative normalisation(s) of the Schechter function (phi_star)}
  \item{kneelo}{a lower bound on the knee parameter}
  \item{slopelo}{a lower bound on the slope parameter}
  \item{normlo}{a lower bound on the norm parameter}
  \item{kneehi}{an upper bound on the knee parameter}
  \item{slopehi}{an upper bound on the slope parameter}
  \item{normhi}{an upper bound on the norm parameter}
  \item{fixk1}{fix the first knee?}
  \item{fixs1}{fix the first slope?}
  \item{fixn1}{fix the first normalisation?}
  \item{fixk2}{fix the second knee?}
  \item{fixs2}{fix the second slope?}
  \item{fixn2}{fix the second normalisation?}
  \item{range}{data range of interest}
  \item{lim1}{lower data limit for fitting}
  \item{lim2}{upper data limit for fitting}
  \item{numlim}{lower number (per bin) limit for fitting}
  \item{method}{choice of 'nlminb' (recommended) or one of 'optim's minimisation methods (e.g., 'Nelder-Mead')}
  \item{volume}{total volume across which the data has been collected (default 1 if vmax = NA)}
  \item{bw}{bin width sizes}
  \item{mag}{are the input data magnitudes?}
  \item{log}{are the input data logged?}
  \item{null}{value of null}
  \item{error}{parameter error estimation method [jack/boot]}
  \item{subvol}{number of sub-volumes to split the input data into}
  \item{sampnum}{number of samplings to be made for bootstrapping method}
  \item{msun}{absolute magnitude of the sun (to be used in calculation of the luminosity density should mag=TRUE)}
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
  \item{j}{luminosity density}
  \item{jlo}{lower luminosity density limit}
  \item{jhi}{upper luminosity density limit}
  \item{chi2}{chi2 value}
  \item{dof}{number of degrees of freedom}
  \item{rchi2}{reduced chi2}
  \item{pval}{probability of observing chi2 value by chance [significant: pval < 0.05]}
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

\keyword{data}
