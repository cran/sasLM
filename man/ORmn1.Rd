\name{ORmn1}
\alias{ORmn1}
\title{Odds Ratio and Score CI of two groups without strata by the MN method}
\description{
  Odds ratio and its score confidence interval of two groups without stratification
}
\usage{
  ORmn1(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
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
  It calculates odds ratio and its score confidence interval of the two groups. The confidence interval is asymmetric, and there is no standard error in the output. This does not support stratification. This implementation uses uniroot function, which usually gives at least 5 significant digits. Whereas PropCIs::orscoreci function uses incremental or decremental search by the factor of 1.001 which gives only less than 3 significant digits.
}
\value{
  There is no standard error.
  \item{odd1}{odd from the first group, y1/(n1 - y1)}
  \item{odd2}{odd from the second group, y2/(n2 - y2)}
  \item{OR}{odds ratio, odd1/odd2}
  \item{lower}{lower confidence limit of OR}
  \item{upper}{upper confidence limit of OR}
}
\examples{
  ORmn1(104, 11037, 189, 11034)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}

\author{Kyun-Seop Bae k@acr.kr}
