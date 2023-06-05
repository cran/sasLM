\name{RRmn}
\alias{RRmn}
\title{Relative Risk and Score CI of two groups with strata by the MN method}
\description{
  Relative risk and its score confidence interval of two groups with stratification by the Miettinen and Nurminen method
}
\usage{
  RRmn(d0, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a strata. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each stratum. The second group is usually the control group.}
  \item{conf.level}{confidence level}
  \item{eps}{absolute value less than eps is regarded as negligible}
}
\details{
  It calculates relative risk and its score confidence interval of the two groups. The confidence interval is asymmetric, and there is no standard error in the output. This supports stratification. This implementation uses uniroot function, which usually gives at least 5 significant digits. Whereas PropCIs::riskscoreci function uses cubic equation approximation which gives only about 2 significant digits. This can be used for meta-analysis also.
}
\value{
  The following output will be returned for each strata and common value. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RR}{relative risk, p1/p2. Point estimate of common RR is calculated with MN weight.}
  \item{lower}{lower confidence limit of RR}
  \item{upper}{upper confidence limit of RR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  RRmn(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{ORmn}}, \code{\link{RDinv}}, \code{\link{RRinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\references{
  Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}
\keyword{2x2 Table}
\author{Kyun-Seop Bae k@acr.kr}
