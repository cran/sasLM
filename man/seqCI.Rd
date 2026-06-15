\name{seqCI}
\alias{seqCI}
\title{Confidence interval with the last Z-value for a group sequential design}
\description{
  Confidence interval with given upper bounds, times of tests, the last Z-value, and confidence level.
}
\usage{
  seqCI(bi, ti, Zval, conf.level=0.95)
}
\arguments{
  \item{bi}{upper bound z-values}
  \item{ti}{times for the tests. These should be in [0, 1].}
  \item{Zval}{the last z-value from the observed data. This is not necessarily the planned final Z-value.}
  \item{conf.level}{confidence level}
}
\details{
  It calculates the confidence interval with given upper bounds, times of tests, the last Z-value, and confidence level. It assumes a two-sided test. \code{mvtnorm::pmvt} (with noncentrality) is better than \code{pmvnorm} in calculating power, sample size, and confidence interval. However, Lan-DeMets used the multivariate normal rather than the multivariate noncentral t distribution. This function follows Lan-DeMets for consistency with previous results. For the theoretical background, see the reference.
}
\value{
  confidence interval of the Z-value for the given confidence level.
}
\examples{
  seqCI(bi = c(2.53, 2.61, 2.57, 2.47, 2.43, 2.38), 
        ti = c(.2292, .3333, .4375, .5833, .7083, .8333), Zval=2.82)
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
