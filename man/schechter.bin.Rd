\name{schechter.bin}
\alias{schechter.bin}

\title{Calculate Schechter Function Binned Variables}

\description{This function calculates a Schechter function fit to a set of input data.}

\usage{schechter.bin(data, vmax = NA, range = range(data), lim1 = NA, 
              lim2 = NA, numlim = 1, volume = max(vmax), bw = 0.1, 
              null = 1E-9)}

\arguments{
  \item{data}{input data vector}
  \item{vmax}{vector of maximum comoving volumes within which object could lie}
  \item{range}{data range of interest}
  \item{lim1}{lower data limit for fitting}
  \item{lim2}{upper data limit for fitting}
  \item{numlim}{lower number (per bin) limit for fitting}
  \item{volume}{total volume across which the data has been collected (default 1 if vmax = NA)}
  \item{bw}{bin width sizes}
  \item{null}{value of null}
}

\value{
  \item{bins}{bin boundaries}
  \item{binmid}{bin midpoints}
  \item{binlo}{bin lower boundaries}
  \item{binhi}{bin upper boundaries}
  \item{num}{number per bin (note: if using vmax values, these are weights)}
  \item{den}{density per bin}
  \item{err}{error per bin}
  \item{errlo}{lower error limit}
  \item{errhi}{upper error limit}
  \item{fitbinmid}{bin midpoints after data limits applied}
  \item{fitbinlo}{bin lower boundaries after data limits applied}
  \item{fitbinhi}{bin upper boundaries after data limits applied}
  \item{fitnum}{number per bin after data limits applied (note: if using vmax values, these are weights)}
  \item{fitden}{density per bin after data limits applied}
  \item{fiterr}{error per bin after data limits applied}
  \item{fiterrlo}{lower error limit after data limits applied}
  \item{fiterrhi}{upper error limit after data limits applied}
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
