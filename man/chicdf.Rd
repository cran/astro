\name{chicdf}
\alias{chicdf}

\title{The Cumulative Distribution Function for the Chi-Squared Distribution}

\description{Calculates the CDF for the chi-squared distribution.}

\usage{chicdf(X, k)}

\arguments{
  \item{X}{a vector of input chi-squared values}
  \item{k}{the number of degrees of freedom}
}

\details{
The chi-squared distribution is the sum of the squares of k independent standard normal random variables, where k represents the number of degrees of freedom. Typically, k is estimated using the relation k = N - n, where N represents the number of data points (observations) in your data set, and n represents the number of fitted parameters in your model.

This function returns the cumulative distribution for a vector of given chi-squared values with an associated number k degrees of freedom.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
X = seq(0,8,len=1000)
kvals = c(1,2,3,4,6,9)
cols = c("yellow","green","turquoise","blue","purple","red")
aplot(NA, type="n", xlab=bquote(chi^2), 
ylab=bquote(paste(F[k],"(",chi^2,")",sep="")), xlim=c(0,8), 
ylim=c(0,1), main="Chi-Squared Cumulative Distribution Function")
grid(lty=1, col="grey90")
for(i in 1:length(kvals)){
    lines(X, chicdf(X=X, k=kvals[i]), lwd=3, col=cols[i])
}
legend("bottomright", legend=paste("k =",kvals), lwd=3, col=cols, 
bty="o", bg=acol("white",alpha=0.7), inset=0.04, text.font=3, 
box.col=NA)
}

\keyword{data}
