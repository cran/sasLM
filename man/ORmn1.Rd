\name{ORmn1}
\alias{ORmn1}
\title{Odds Ratio and Score CI of two groups}
\description{
  Odds ratio and its score confidence interval of two groups without stratification
}
\usage{
  ORmn1(y1, n1, y2, n2, conf.level=0.95)
}
\arguments{
  \item{y1}{positive event count of test (the first) group}
  \item{n1}{total count of the test (the first) group}
  \item{y2}{positive event count of control (the second) group}
  \item{n2}{total count of control (the second) group}
  \item{conf.level}{confidence level}
}
\details{
  It calculates odds ratio and its score confidence interval of two groups. The confidence interval is asymmetric and there is no standard error in the output. This does not support stratification. This implementation uses nlminb minimization function.
}
\value{
  The result is a data.frame. There is no standard error.
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
  \code{\link{OR}}, \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\references{
Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26
}

\author{Kyun-Seop Bae k@acr.kr}
