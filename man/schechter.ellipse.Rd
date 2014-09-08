\name{schechter.ellipse}
\alias{schechter.ellipse}
\alias{.schechter.ellipse.fit}

\title{Calculate Error Matrices for a Schechter Function Fit}

\description{This function calculates ellipsoidal error matrices for a given Schechter Function fit.}

\usage{schechter.ellipse(data, vmax = NA, knee, slope, norm, chi2, 
                  datarange = NA, kneerange = c(-24,-16), 
                  sloperange = c(-2,1.5), kneeofflims = NA, 
                  slopeofflims = NA, kneestep = 0.5, slopestep = 0.1, 
                  kneesteps = NA, slopesteps = NA, lim1 = NA, lim2 = NA, 
                  numlim = 1, method = "nlminb", volume = max(vmax), 
                  bw = 0.1, mag = FALSE, log = FALSE, null = 1E-9, ...)}

\arguments{
  \item{data}{input data vector}
  \item{vmax}{vector of maximum comoving volumes within which object could lie}
  \item{knee}{the knees(s) of the Schechter function (L_star/M_star)}
  \item{slope}{the slope(s) of the Schechter function (alpha)}
  \item{norm}{the normalisation(s) of the Schechter function (phi_star)}
  \item{chi2}{the full chi2 result from this fit (not reduced)}
  \item{datarange}{the range across which the data is evaluated}
  \item{kneerange}{range of knee values}
  \item{sloperange}{range of slope values}
  \item{kneeofflims}{alternative to kneerange, vector length 2 describing limit offsets from knees}
  \item{slopeofflims}{alternative to sloperange, vector length 2 describing limit offsets from slopes}
  \item{kneestep}{the matrix step in knee values}
  \item{slopestep}{the matrix step in slope values}
  \item{kneesteps}{alternative to kneestep, the number of steps in the matrix}
  \item{slopesteps}{alternative to slopestep, the number of steps in the matrix}
  \item{lim1}{lower data limit for fitting}
  \item{lim2}{upper data limit for fitting}
  \item{numlim}{lower number (per bin) limit for fitting}
  \item{method}{choice of 'nlminb' (recommended) or one of 'optim's minimisation methods (e.g., 'Nelder-Mead')}
  \item{volume}{total volume across which the data has been collected (default 1 if vmax = NA)}
  \item{bw}{bin width sizes}
  \item{mag}{are the input data magnitudes?}
  \item{log}{are the input data logged?}
  \item{null}{value of null}
  \item{...}{additional arguments to be passed to 'integrate'}
}

\value{
  \item{knees}{knee bin midpoints}
  \item{slopes}{slope bin midpoints}
  \item{res1}{result matrix for knee1/slope1}
  \item{res2}{result matrix for knee2/slope2}
  \item{s1}{1 sigma chi2 limit (min+2.30)}
  \item{s2}{2 sigma chi2 limit (min+6.17)}
  \item{s3}{3 sigma chi2 limit (min+11.8)}
  \item{kneelo1}{lower knee error (all sigmas)}
  \item{kneehi1}{upper knee error (all sigmas)}
  \item{slopelo1}{lower slope error (all sigmas)}
  \item{slopehi1}{upper slope error (all sigmas)}
  \item{kneelo2}{lower secondary knee error (all sigmas)}
  \item{kneehi2}{upper secondary knee error (all sigmas)}
  \item{slopelo2}{lower secondary slope error (all sigmas)}
  \item{slopehi2}{upper secondary slope error (all sigmas)}
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
