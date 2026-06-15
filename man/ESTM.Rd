\name{ESTM}
\alias{ESTM}
\alias{ESTIMATE}
\title{Estimate Linear Function}
\description{
  Estimates Linear Function with a formula and a dataset.
}
\usage{
  ESTM(L, Formula, Data, conf.level=0.95)
}
\arguments{
  \item{L}{a matrix of linear function rows to be tested}
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{conf.level}{confidence level of the confidence limit}
}
\details{
  It tests rows of linear functions. A linear function means a linear combination of estimated coefficients. It is similar to the ESTIMATE statement of SAS PROC GLM. This is a convenient version of the \code{est} function.
}
\value{
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for the t distribution}
  \item{Df}{degrees of freedom}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual degrees of freedom}
}
\examples{
  ESTM(t(c(0, -1, 1)), uptake ~ Type, CO2) # Quevec - Mississippi 
}
\seealso{\code{\link{est}}}
\author{Kyun-Seop Bae k@acr.kr}
