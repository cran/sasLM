\name{RRmn}
\alias{RRmn}
\title{Relative Risk and Score CI of two groups with strata}
\description{
  Relative risk and its score confidence interval of two groups with stratification
}
\usage{
  RRmn(d0, conf.level=0.95, weight="MN", eps=1e-8)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a strata. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each strata. The second group is usually the control group.}
  \item{conf.level}{confidence level}
  \item{weight}{Default weight is "MN"(Miettinen-Nurminen) method. If not, use "MH"(Mantel-Haenszel) method.}
  \item{eps}{absolute value less than eps is regarded as 0}
}
\details{
  It calculates relative risk and its score confidence interval of two groups. The confidence interval is asymmetric and there is no standard error in the output. This supports stratification. This implementation uses nlminb minimization function.
}
\value{
  The following output will be returned for each strata and common value. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RR}{relative risk, p1/p2. Point estimate of common RR is calculated with MH weight.}
  \item{lower}{lower confidence limit of RR}
  \item{upper}{upper confidence limit of RR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  RRmn(d1)
}
\seealso{
  \code{\link{RR}}, \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{ORmn}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}
\author{Kyun-Seop Bae k@acr.kr}
