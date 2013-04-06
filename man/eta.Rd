\name{eta}
\alias{eta}

\title{Calculate ETA of a For-Loop}

\description{Calculates an ETA of a for-loop completing the loop, based on how many loops are left, and the start time of the first loop.}

\usage{eta(i, total, start)}

\arguments{
  \item{i}{current element}
  \item{total}{total number of elements}
  \item{start}{start processor time (seconds)}
}

\details{
The start time is best given by 'proc.time()[3]'.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
