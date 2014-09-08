\name{p2chi}
\alias{p2chi}

\title{The Chi-Squared Statistic for a Given P-Value}

\description{Calculates the chi-squared statistic for a given p-value.}

\usage{p2chi(p, k, steps = 30, chimax = 1E31)}

\arguments{
  \item{p}{a vector of input p-values}
  \item{k}{the number of degrees of freedom}
  \item{steps}{number of accuracy steps}
  \item{chimax}{the upper chi-squared search limit}
}

\details{
The chi-squared distribution is the sum of the squares of k independent standard normal random variables, where k represents the number of degrees of freedom. Typically, k is estimated using the relation k = N - n, where N represents the number of data points (observations) in your data set, and n represents the number of fitted parameters in your model.

The p-value represents the estimated probability of rejecting the null hypothesis. Here we assume the null hypothesis to be that the sample follows a chi-squared distribution as parameterised by the number of degrees of freedom k. Typically, if p >= 0.05, the data appear to be consistent with the null hypothesis, and if p < 0.05, there is significant evidence against the null hypothesis in favour of the alternative. These limits, albeit rather arbitrary, are nevertheless consistently used in the literature. Often, the p-value is (incorrectly) interpreted as the probability that the null hypothesis is true.

This function returns the chi-squared statistic for a vector of given p-values with an associated number k degrees of freedom. Increasing the number of accuracy steps will increase the final accuracy, however; steps = 30 (the default) should be sufficient for most tasks.
}

\author{Lee Kelvin <lee.kelvin@uibk.ac.at>}

\seealso{

The astronomy package: \code{\link{astro}}.

}

\examples{
X = seq(0,20,len=1000)
kvals = c(1,2,3,4,6,9)
cols = c("yellow","green","turquoise","blue","purple","red")
for(i in 1:length(kvals)){
    psig = p2chi(p=0.05, k=kvals[i])
    pgood = seq(0,psig,len=100)
    pbad = seq(psig,max(X),len=100)
    aplot(X, chipdf(X=X, k=kvals[i]), type="n", xlab=bquote(chi^2), 
    ylab=bquote(paste(f[k],"(",chi^2,")",sep="")), 
    ylim=c(0,min(0.5,max(chipdf(X=X, k=kvals[i])))), 
    main="Chi-Squared Probability Density Function")
    polygon(x=c(pgood,rev(pgood)), 
    y=c(chipdf(pgood, k=kvals[i]),rep(0,len=length(pgood))), 
    col=acol("yellow",alpha=0.5), border=NA)
    polygon(x=c(pbad,rev(pbad)), 
    y=c(chipdf(pbad, k=kvals[i]),rep(0,len=length(pbad))), 
    col=acol("turquoise",alpha=0.5), border=NA)
    lines(X, chipdf(X=X, k=kvals[i]), lwd=3, col=cols[i])
    legend("topright", 
    legend=c("p >= 0.05 (data consistent with null hypothesis)", 
    "p < 0.05 (reject null hypothesis)"), bty="o", 
    bg=acol("white",alpha=0.7), inset=0.04, text.font=3, box.col=NA, 
    fill=acol(c("yellow","turquoise"),alpha=0.5), border=NA)
    legend("right", legend=paste("k =",kvals[i]), lwd=3, col=cols[i], 
    bty="o", bg=acol("white",alpha=0.7), inset=0.04, text.font=3, 
    box.col=NA)
}
}

\keyword{data}
