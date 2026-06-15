\name{ScoreCI}
\alias{ScoreCI}
\title{Score Confidence Interval for a Proportion or a Binomial Distribution}
\description{
  Score confidence interval of a proportion in one group
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
  It calculates the score confidence interval of a proportion in one group. The confidence interval is asymmetric, and there is no standard error in the output. If you need the percent scale, multiply the output by 100.
}
\value{
  The result is a \code{data.frame}. There is no standard error.
  \item{PE}{point estimate of the proportion}
  \item{Lower}{lower confidence limit of the proportion}
  \item{Upper}{upper confidence limit of the proportion}
}
\examples{
  ScoreCI(104, 11037)
}
\seealso{
  \code{\link{binom.test}}, \code{\link{prop.test}}
}
\author{Kyun-Seop Bae k@acr.kr}
