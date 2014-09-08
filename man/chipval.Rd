\name{chipval}
\alias{chipval}

\title{The P-Value for the Chi-Squared Distribution}

\description{Calculates the p-value for a given chi-squared statistic.}

\usage{chipval(X, k)}

\arguments{
  \item{X}{a vector of input chi-squared values}
  \item{k}{the number of degrees of freedom}
}

\details{
The chi-squared distribution is the sum of the squares of k independent standard normal random variables, where k represents the number of degrees of freedom. Typically, k is estimated using the relation k = N - n, where N represents the number of data points (observations) in your data set, and n represents the number of fitted parameters in your model.

The p-value represents the estimated probability of rejecting the null hypothesis. Here we assume the null hypothesis to be that the sample follows a chi-squared distribution as parameterised by the number of degrees of freedom k. Typically, if p >= 0.05, the data appear to be consistent with the null hypothesis, and if p < 0.05, there is significant evidence against the null hypothesis in favour of the alternative. These limits, albeit rather arbitrary, are nevertheless consistently used in the literature. Often, the p-value is (incorrectly) interpreted as the probability that the null hypothesis is true.

This function returns the p-value for a vector of given chi-squared values with an associated number k degrees of freedom.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
X = seq(0,20,by=1)
kvals = c(1,2,3,4,6,9)
par("mar"=c(1,1,3,1))
plot(NA, type="n", xlim=c(-0.5,6.5), ylim=c(21.5,-0.5), axes=FALSE, 
xlab="", ylab="", xaxs="i", yaxs="i", main="P-Value Lookup Table")
text(x=0, y=0, labels=bquote(chi^2))
for(i in 1:length(kvals)){
    text(x=i, y=0, labels=paste("k =",kvals[i]))
}
abline(v=seq(0.5,length(kvals)-0.5,by=1), col="grey90", lwd=3)
for(i in 1:length(X)){
    text(x=0, y=i, labels=paste(X[i]))
}
abline(h=seq(0.5,length(X)-0.5,by=1), col="grey90", lwd=3)
cols = acol(colorRampPalette(c("green","red"))(100),alpha=0.5)
for(i in 1:length(X)){
    for(j in 1:length(kvals)){
        p = chipval(X=X[i], k=kvals[j])
        col = ((p^0.35)*(length(cols)-1))+1
        rect(xleft=j-0.5, xright=j+0.5, ybottom=i-0.5, ytop=i+0.5, 
        col=cols[col], border=NA)
        text(x=j, y=i, labels=formatC(p,format="f",digits=3))
    }
}
par("mar"=c(5.1,4.1,4.1,2.1))
}

\keyword{data}
