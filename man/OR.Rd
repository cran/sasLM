\name{OR}
\alias{OR}
\title{Odds Ratio of two groups}
\description{
  Odds ratio between two groups
}
\usage{
  OR(y1, n1, y2, n2, conf.level=0.95)
}
\arguments{
  \item{y1}{positive event count of the test (the first) group}
  \item{n1}{total count of the test (the first) group}
  \item{y2}{positive event count of the control (the second) group}
  \item{n2}{total count of the control (the second) group}
  \item{conf.level}{confidence level}
}
\details{
  It calculates the odds ratio of two groups. No continuity correction is done here. If you need the percent scale, multiply the output by 100.
}
\value{
  The result is a \code{data.frame}.
  \item{odd1}{odds from the first group, y1/(n1 - y1)}
  \item{odd2}{odds from the second group, y2/(n2 - y2)}
  \item{OR}{odds ratio, odd1/odd2}
  \item{SElog}{standard error of log(OR)}
  \item{lower}{lower confidence limit of OR}
  \item{upper}{upper confidence limit of OR}
}
\examples{
  OR(104, 11037, 189, 11034) # no continuity correction
}
\seealso{
  \code{\link{RD}}, \code{\link{RR}}, \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\author{Kyun-Seop Bae k@acr.kr}
