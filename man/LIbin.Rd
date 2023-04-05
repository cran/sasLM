\name{LIbin}
\alias{LIbin}
\title{Likelihood Interval for a Proportion or a Binomial Distribution}
\description{
  Likelihood interval of a proportion in one group
}
\usage{
  LIbin(y, n, k, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{y}{positive event count of a group}
  \item{n}{total count of a group}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
  \item{eps}{Values less than eps are considered as 0.}
}
\details{
  It calculates likelihood interval of a proportion in one group. The likelihood interval is asymmetric and there is no standard error in the output. If you need percent scale, multiply the output by 100.
}
\value{
  \item{y}{positive (concerning) event count}
  \item{n}{total trial count}
  \item{PE}{point estimation for the proportion}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  LIbin(3, 14, k=2)
  LIbin(3, 14, k=5)
  LIbin(3, 14, k=15)
  LIbin(3, 14)
  ScoreCI(3, 14)
  binom.test(3, 14)
  prop.test(3, 14)
}
\seealso{
  \code{\link{binom.test}}, \code{\link{prop.test}}, \code{\link{ScoreCI}}
}
\references{
  Fisher RA. Statistical methods and scientific inference. 3e. 1973. pp68-76.  
}
\author{Kyun-Seop Bae k@acr.kr}
