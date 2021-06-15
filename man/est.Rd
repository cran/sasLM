\name{est}
\alias{est}
\title{Estimate Linear Functions}
\description{
  Estimates Linear Functions with a given GLM result.
}
\usage{
  est(L, X, rx, conf.level=0.95, adj="lsd")
}
\arguments{
  \item{L}{a matrix of linear contrast rows to be tested}
  \item{X}{a model (design) matrix from \code{ModelMatrix}}
  \item{rx}{a result of \code{lfit} function}
  \item{conf.level}{confidence level of confidence limit}
  \item{adj}{"lsd" or "tukey" to adjust p-value and confidence limit}
}
\details{
  It tests rows of linear function. Linear function means linear combination of estimated coefficients. It corresponds to SAS PROC GLM ESTIMATE. Same sample size per group is assumed for the Tukey adjustment.
}
\value{
  \item{Estimate}{point estimate of the input linear constrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Df}{degree of freedom}
  \item{Pr(>|t|)}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
  x = ModelMatrix(uptake ~ Type, CO2)
  rx = REG(uptake ~ Type, CO2, summarize=FALSE)
  est(t(c(0, -1, 1)), x$X, rx) # Quevec - Mississippi 
  t.test(uptake ~ Type, CO2) # compare with the above
}
\seealso{\code{\link{ESTM}}}
\author{Kyun-Seop Bae k@acr.kr}
