\name{ScoreCI}
\alias{ScoreCI}
\title{Score CI of a proportion by Wilson}
\description{
  Score confidence of a proportion in one group by Wilson
}
\usage{
  ScoreCI(y, n, conf.level=0.95)
}
\arguments{
  \item{y}{positive event count of a group}
  \item{n}{total count of a group}
  \item{conf.level}{confidence level}
}
\details{
  It calculates score confidence interval of a proportion in one group as Wilson suggested. If you need percent scale, multiply the output by 100.
}
\value{
  The result is a data.frame. There is no standard error.
  \item{Prop}{proportion from the first group}
  \item{lower}{lower confidence limit of Prop}
  \item{upper}{upper confidence limit of Prop}
}
\examples{
  ScoreCI(104, 11037)
}
\seealso{
  \code{\link{binom.test}}
}
\author{Kyun-Seop Bae k@acr.kr}
