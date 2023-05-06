\name{ScoreCI}
\alias{ScoreCI}
\title{Score Confidence Interval for a Proportion or a Binomial Distribution}
\description{
  Score confidence of a proportion in one group
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
  It calculates score confidence interval of a proportion in one group. The confidence interval is asymmetric and there is no standard error in the output. If you need percent scale, multiply the output by 100.
}
\value{
  The result is a data.frame. There is no standard error.
  \item{PE}{point estimation for the proportion}
  \item{Lower}{lower confidence limit of Prop}
  \item{Upper}{upper confidence limit of Prop}
}
\examples{
  ScoreCI(104, 11037)
}
\seealso{
  \code{\link{binom.test}}, \code{\link{prop.test}}
}
\author{Kyun-Seop Bae k@acr.kr}
