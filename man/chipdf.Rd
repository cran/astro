\name{chipdf}
\alias{chipdf}

\title{The Probability Density Function for the Chi-Squared Distribution}

\description{Calculates the PDF for the chi-squared distribution.}

\usage{chipdf(X, k)}

\arguments{
  \item{X}{a vector of input chi-squared values}
  \item{k}{the number of degrees of freedom}
}

\details{
The chi-squared distribution is the sum of the squares of k independent standard normal random variables, where k represents the number of degrees of freedom. Typically, k is estimated using the relation k = N - n, where N represents the number of data points (observations) in your data set, and n represents the number of fitted parameters in your model.

This function returns the probability density function for a vector of given chi-squared values with an associated number k degrees of freedom.
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
ylab=bquote(paste(f[k],"(",chi^2,")",sep="")), xlim=c(0,8), 
ylim=c(0,0.5), main="Chi-Squared Probability Density Function")
grid(lty=1, col="grey90")
for(i in 1:length(kvals)){
    lines(X, chipdf(X=X, k=kvals[i]), lwd=3, col=cols[i])
}
legend("topright", legend=paste("k =",kvals), lwd=3, col=cols, 
bty="o", bg=acol("white",alpha=0.7), inset=0.04, text.font=3, 
box.col=NA)
}

\keyword{data}
