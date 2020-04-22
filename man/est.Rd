\name{est}
\alias{est}
\title{Estimate Linear Contrast}
\description{
  Estimates Linear Contrast(s) with a given GLM result.
}
\usage{
  est(L, rx)
}
\arguments{
  \item{L}{a matrix of linear contrast rows to be tested}
  \item{rx}{a result of \code{lfit} function}
}
\details{
  It tests rows of linear contrast. It corresponds to SAS PROC GLM ESTIMATE.
}
\value{
  Returns a table of expectations, t values and p-values.
  \item{Estimate}{point estimate of the input linear constrast}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Pr(>|t|}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
  x = ModelMatrix(uptake ~ Type, CO2)
  y = model.frame(uptake ~ Type, CO2)[,1]
  rx = lfit(x, y)
  est(t(c(0, -1, 1)), rx) # Quevec - Mississippi 
  t.test(uptake ~ Type, CO2) # compare with the above
}
\author{Kyun-Seop Bae k@acr.kr}
