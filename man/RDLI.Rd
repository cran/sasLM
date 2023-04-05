\name{RDLI}
\alias{RDLI}
\title{Risk (Proportion) Difference and its Likelihood Interval between two groups without strata}
\description{
  Risk difference and its likelihood interval between two groups without stratification
}
\usage{
  RDLI(y1, n1, y2, n2, conf.level=0.95, k, eps=1e-8)
}
\arguments{
  \item{y1}{positive event count of test (the first) group}
  \item{n1}{total count of the test (the first) group. Maximum allowable value is 1e8.}
  \item{y2}{positive event count of control (the second) group}
  \item{n2}{total count of control (the second) group. Maximum allowable value is 1e8.}
  \item{conf.level}{approximate confidence level to calculate k when k is missing.}
  \item{k}{1/k likelihood interval will be provided}
  \item{eps}{absolute value less than eps is regarded as negligible}
}
\details{
  It calculates risk (proportion) difference and its likelihood interval between the two groups. The likelihood interval is asymmetric, and there is no standard error in the output. This does not support stratification.
}
\value{
  There is no standard error.
  \item{p1}{proportion from the first group, y1/n1}
  \item{p2}{proportion from the second group, y2/n2}
  \item{RD}{risk difference, p1 - p2}
  \item{lower}{lower likelihood limit of RD}
  \item{upper}{upper likelihood limit of RD}
}
\examples{
  RDLI(7, 10, 3, 10)
  RD(7, 10, 3, 10)
  RDLI(3, 10, 7, 10)
  RD(3, 10, 7, 10)
}
\seealso{
  \code{\link{RD}}, \code{\link{RDmn1}}
}
\references{
  None yet!
}
\author{Kyun-Seop Bae k@acr.kr}
