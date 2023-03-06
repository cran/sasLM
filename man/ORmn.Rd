\name{ORmn}
\alias{ORmn}
\title{Odds Ratio and Score CI of two groups with strata by MN method}
\description{
  Odds ratio and its score confidence interval of two groups with stratification by the Miettinen and Nurminen method
}
\usage{
  ORmn(d0, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a strata. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each strata. The second group is usually the control group.}
  \item{conf.level}{confidence level}
  \item{eps}{absolute value less than eps is regarded as negligible}
}
\details{
  It calculates odds ratio and its score confidence interval of the two groups. The confidence interval is asymmetric, and there is no standard error in the output. This supports stratification. This implementation uses nlminb minimization function, which usually gives about 7 significant digits. PropCIs::orscoreci function uses incremental or decremental search by the factor of 1.001 which gives only about 3 significant digits.
}
\value{
  The following output will be returned for each stratum and common value. There is no standard error.
  \item{odd1}{odd from the first group, y1/(n1 - y1)}
  \item{odd2}{odd from the second group, y2/(n2 - y2)}
  \item{OR}{odds ratio, odd1/odd2. The point estimate of common OR is calculated with MN weight.}
  \item{lower}{lower confidence limit of OR}
  \item{upper}{upper confidence limit of OR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  ORmn(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{RDinv}}, \code{\link{RRinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}

\author{Kyun-Seop Bae k@acr.kr}
