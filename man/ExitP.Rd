\name{ExitP}
\alias{ExitP}
\title{Exit Probability with Cumulative Z-test in Group Sequential Design}
\description{
  Exit probabilities with the given drift, upper bounds, and times of test.
}
\usage{
  ExitP(Theta, bi, ti=NULL)
}
\arguments{
  \item{Theta}{drift value defined by Lan-DeMets. See the reference.}
  \item{bi}{upper bound z-values}
  \item{ti}{times of test. These should be in the range of [0, 1]. If omitted, even intervals are assumed.}
}
\details{
  It calculates exit probabilities and cumulative exit probabilities with the given drift, upper z-bounds, and times of test. If the times of test are not given, even intervals are assumed. \code{mvtnorm::pmvt} (with noncentrality) is better than \code{pmvnorm} in calculating power and sample size. However, Lan-DeMets used the multivariate normal rather than the multivariate noncentral t distribution. This function follows Lan-DeMets for consistency with previous results.
}
\value{
  The result is a matrix.
  \item{ti}{time of test}
  \item{bi}{upper z-bound}
  \item{cum.alpha}{cumulative alpha-value}
}
\examples{
  b0 = seqBound(ti=(1:5)/5)[, "up.bound"]
  ExitP(Theta = Drift(b0), bi = b0)
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
