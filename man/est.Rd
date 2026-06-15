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
  \item{rx}{a result of the \code{lfit} function}
  \item{conf.level}{confidence level of the confidence limit}
  \item{adj}{adjustment method for grouping. This supports "tukey", "bon", "scheffe", "duncan", and "dunnett". This only affects grouping, not the confidence interval.}
  \item{paired}{If this is \code{TRUE}, the \code{L} matrix is for the pairwise comparison such as that of the \code{PDIFF} function.}
}
\details{
  It tests rows of linear functions. A linear function means a linear combination of estimated coefficients. It corresponds to the ESTIMATE statement of SAS PROC GLM. The same sample size per group is assumed for the Tukey adjustment.
}
\value{
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit by the "lsd" method}
  \item{Upper CL}{upper confidence limit by the "lsd" method}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for the t distribution, for methods other than "scheffe"}
  \item{F value}{value for the F distribution, for the "scheffe" method only}
  \item{Df}{degrees of freedom of the residuals}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual degrees of freedom, for methods other than "scheffe"}
  \item{Pr(>F)}{probability of a larger F value from the F distribution with the residual degrees of freedom, for the "scheffe" method only}
}
\examples{
  x = ModelMatrix(uptake ~ Type, CO2)
  rx = REG(uptake ~ Type, CO2, summarize=FALSE)
  est(t(c(0, -1, 1)), x$X, rx) # Quebec - Mississippi 
  t.test(uptake ~ Type, CO2) # compare with the above
}
\seealso{\code{\link{ESTM}}, \code{\link{PDIFF}}}
\author{Kyun-Seop Bae k@acr.kr}
