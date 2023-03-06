\name{RD}
\alias{RD}
\title{Risk Difference between two groups}
\description{
  Risk (proportion) difference between two groups
}
\usage{
  RD(y1, n1, y2, n2, conf.level=0.95)
}
\arguments{
  \item{y1}{positive event count of test (the first) group}
  \item{n1}{total count of the test (the first) group}
  \item{y2}{positive event count of control (the second) group}
  \item{n2}{total count of control (the second) group}
  \item{conf.level}{confidence level}
}
\details{
  It calculates risk difference between the two groups. No continuity correction here. If you need percent scale, multiply the output by 100.
}
\value{
  The result is a data.frame.
  \item{p1}{proportion from the first group}
  \item{p2}{proportion from the second group}
  \item{RD}{risk difference, p1 - p2}
  \item{SE}{standard error of RD}
  \item{lower}{lower confidence limit of RD}
  \item{upper}{upper confidence limit of RD}
}
\examples{
  RD(104, 11037, 189, 11034) # no continuity correction
}
\seealso{
  \code{\link{RR}}, \code{\link{OR}}, \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}
}
\author{Kyun-Seop Bae k@acr.kr}
