\name{RRinv}
\alias{RRinv}
\title{Relative Risk of two groups with strata by inverse variance method}
\description{
  Relative risk and its score confidence interval of two groups with stratification by inverse variance method
}
\usage{
  RRinv(d0, conf.level=0.95)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a stratum. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each stratum. The second group is usually the control group.}
  \item{conf.level}{confidence level}
}
\details{
  It calculates relative risk and its confidence interval of two groups by inverse variance method. This supports stratification. This can be used for meta-analysis also.
}
\value{
  The following output will be returned for each stratum and common value. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RR}{relative risk, p1/p2. The point estimate of common RR is calculated with MH weight.}
  \item{lower}{lower confidence limit of RR}
  \item{upper}{upper confidence limit of RR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  RRinv(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}, \code{\link{RDinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\author{Kyun-Seop Bae k@acr.kr}
