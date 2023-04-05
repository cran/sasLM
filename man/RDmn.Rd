\name{RDmn}
\alias{RDmn}
\title{Risk Difference and Score CI between two groups with strata by the MN method}
\description{
  Risk difference and its score confidence interval between two groups with stratification by the Miettinen and Nurminen method
}
\usage{
  RDmn(d0, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a stratum. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each stratum. The second group is usually the control group. Maximum allowable value for n1 and n2 is 1e8.}
  \item{conf.level}{confidence level}
  \item{eps}{absolute value less than eps is regarded as negligible}
}
\details{
  It calculates risk difference and its score confidence interval between the two groups. The confidence interval is asymmetric, and there is no standard error in the output. If you need percent scale, multiply the output by 100. This supports stratification. This implementation uses uniroot function which usually gives at least 5 significant digits. This can be used for meta-analysis also.
}
\value{
  The following output will be returned for each stratum and common value. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RD}{risk difference, p1 - p2. The point estimate of common RD is calculated with MN weight.}
  \item{lower}{lower confidence limit of RD}
  \item{upper}{upper confidence limit of RD}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  RDmn(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RRmn}}, \code{\link{ORmn}}, \code{\link{RDinv}}, \code{\link{RRinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\references{
  Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}
\author{Kyun-Seop Bae k@acr.kr}
