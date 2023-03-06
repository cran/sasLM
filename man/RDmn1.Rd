\name{RDmn1}
\alias{RDmn1}
\title{Risk Difference and Score CI between two groups without strata by the MN method}
\description{
  Risk difference and its score confidence interval between two groups without stratification
}
\usage{
  RDmn1(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
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
  It calculates risk difference and its score confidence interval between the two groups. The confidence interval is asymmetric, and there is no standard error in the output. If you need percent scale, multiply the output by 100. This does not support stratification. This implementation uses nlminb minimization function.
}
\value{
  The result is a data.frame. There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RD}{risk difference, p1 - p2}
  \item{lower}{lower confidence limit of RD}
  \item{upper}{upper confidence limit of RD}
}
\examples{
  RDmn1(104, 11037, 189, 11034)
}
\seealso{
  \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}
\author{Kyun-Seop Bae k@acr.kr}
