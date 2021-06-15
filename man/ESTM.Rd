\name{ESTM}
\alias{ESTM}
\alias{ESTMIMATE}
\title{Estimate Linear Function}
\description{
  Estimates Linear Function with a formula and a dataset.
}
\usage{
  ESTM(L, Formula, Data, conf.level=0.95)
}
\arguments{
  \item{L}{a matrix of linear functions rows to be tested}
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{conf.level}{confidence level of confidence limit}
}
\details{
  It tests rows of linear functions. Linear function means linear combination of estimated coefficients. It is similar to SAS PROC GLM ESTIMATE. This is a little convenient version of \code{est} function.
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
  ESTM(t(c(0, -1, 1)), uptake ~ Type, CO2) # Quevec - Mississippi 
}
\seealso{\code{\link{est}}}
\author{Kyun-Seop Bae k@acr.kr}
