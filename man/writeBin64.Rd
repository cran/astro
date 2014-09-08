\name{writeBin64}
\alias{writeBin64}

\title{Transfer Binary Data To a Connection}

\description{Write binary data to a connection.}

\usage{writeBin64(object, con, size = NA_integer_,
           endian = .Platform$endian, useBytes = FALSE)}

\arguments{
  \item{object}{An R object to be written to the connection.}
  \item{con}{A connection object or a character string naming a file or a raw vector.}
  \item{size}{integer.  The number of bytes per element in the byte stream. The default, 'NA_integer_', uses the natural size. Size changing is not supported for raw and complex vectors.}
  \item{endian}{The endian-ness ('big' or 'little') of the target system for the file.  Using 'swap' will force swapping endian-ness.}
  \item{useBytes}{See 'writeLines'.}
}

\details{
An extension of the standard R function \code{\link{writeBin}}. This 64-bit extension splits data into 2^30 chunks, and writes each out sequentially, overcoming the 2^31-1 limit of the original function.
}

\author{
Aaron Robotham

Maintainer: Lee Kelvin <lee.kelvin@uibk.ac.at>
}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\keyword{data}
