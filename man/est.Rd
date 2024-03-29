\name{est}
\alias{est}
\title{Estimate Linear Functions}
\description{
  Estimates Linear Functions with a given GLM result.
}
\usage{
  est(L, X, rx, conf.level=0.95, adj="lsd", paired=FALSE)
}
\arguments{
  \item{L}{a matrix of linear contrast rows to be tested}
  \item{X}{a model (design) matrix from \code{ModelMatrix}}
  \item{rx}{a result of \code{lfit} function}
  \item{conf.level}{confidence level of confidence limit}
  \item{adj}{adjustment method for grouping. This supports "tukey", "bon", "scheffe", "duncan", and "dunnett". This only affects grouping, not the confidence interval.}
  \item{paired}{If this is TRUE, L matrix is for the pairwise comparison such as PDIFF function.}
}
\details{
  It tests rows of linear function. Linear function means linear combination of estimated coefficients. It corresponds to SAS PROC GLM ESTIMATE. Same sample size per group is assumed for the Tukey adjustment. 
}
\value{
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit by "lsd" method}
  \item{Upper CL}{upper confidence limit by "lsd" method}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution for other than "scheffe" method}
  \item{F value}{value for F distribution for "scheffe" method only}
  \item{Df}{degree of freedom of residuals}
  \item{Pr(>|t|)}{probability of larger than absolute t value from t distribution with residual's degree of freedom, for other than "scheffe" method}
  \item{Pr(>F)}{probability of larger than F value from F distribution with residual's degree of freedom, for "scheffe" method only}
}
\examples{
  x = ModelMatrix(uptake ~ Type, CO2)
  rx = REG(uptake ~ Type, CO2, summarize=FALSE)
  est(t(c(0, -1, 1)), x$X, rx) # Quebec - Mississippi 
  t.test(uptake ~ Type, CO2) # compare with the above
}
\seealso{\code{\link{ESTM}}, \code{\link{PDIFF}}}
\author{Kyun-Seop Bae k@acr.kr}
