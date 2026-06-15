\name{RDinv}
\alias{RDinv}
\title{Risk Difference between two groups with strata by inverse variance method}
\description{
  Risk difference and its confidence interval between two groups with stratification by the inverse variance method
}
\usage{
  RDinv(d0, conf.level=0.95)
}
\arguments{
  \item{d0}{A \code{data.frame} or matrix, of which each row represents a stratum. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for the sample size of each stratum. The second group is usually the control group.}
  \item{conf.level}{confidence level}
}
\details{
  It calculates the risk difference and its confidence interval between two groups by the inverse variance method. The common risk difference is given by both the fixed effect model and the DerSimonian-Laird random effects model, with Cochran's Q test for heterogeneity. If you need the percent scale, multiply the output by 100. This can be used for meta-analysis also.
}
\value{
  \item{RDs}{risk difference and its confidence interval of each stratum}
  \item{Heterogeneity}{Cochran's Q statistic for heterogeneity across the strata and its p-value}
  \item{tau2}{between-strata variance estimated by the method of moments}
  \item{Fixed}{common risk difference, its standard error, and its confidence interval by the fixed effect model with the inverse variance weights}
  \item{Random}{common risk difference, its standard error, and its confidence interval by the DerSimonian-Laird random effects model}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  RDinv(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}, \code{\link{RRinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\keyword{Meta-analysis}
\author{Kyun-Seop Bae k@acr.kr}
