\name{RRmn1}
\alias{RRmn1}
\title{Relative Risk and Score CI of two groups without strata by by MN method}
\description{
  Relative risk and its score confidence interval of the two groups without stratification
}
\usage{
  RRmn1(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{y1}{positive event count of test (the first) group}
  \item{n1}{total count of the test (the first) group}
  \item{y2}{positive event count of control (the second) group}
  \item{n2}{total count of control (the second) group}
  \item{conf.level}{confidence level}
  \item{eps}{absolute value less than eps is regarded as negligible}    
}
\details{
  It calculates the relative risk and its score confidence interval of the two groups. The confidence interval is asymmetric, and there is no standard error in the output. This does not support stratification. This implementation uses nlminb minimization function, which usually gives about 7 significant digits. PropCIs::riskscoreci function uses cubic equation approximation which gives only about 2 significant digits.
}
\value{
  The result is a data.frame. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RR}{relative risk, p1/p2}
  \item{lower}{lower confidence limit of RR}
  \item{upper}{upper confidence limit of RR}
}
\examples{
  RRmn1(104, 11037, 189, 11034)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}
\author{Kyun-Seop Bae k@acr.kr}
