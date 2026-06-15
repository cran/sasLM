\name{Drift}
\alias{Drift}
\title{Drift defined by Lan and DeMets for Group Sequential Design}
\description{
  Calculate the drift value with given upper bounds (z-values), times of test, and power.
}
\usage{
  Drift(bi, ti=NULL, Power=0.9)
}
\arguments{
  \item{bi}{upper bound z-values}
  \item{ti}{times of test. These should be in the range of [0, 1]. If omitted, equal intervals are assumed.}
  \item{Power}{target power at the final test}
}
\details{
  It calculates the drift value with given upper bound z-values, times of test, and power. If the times of test are not given, equal intervals are assumed. \code{mvtnorm::pmvt} (with noncentrality) is better than \code{pmvnorm} in calculating power and sample size. However, Lan-DeMets used the multivariate normal rather than the multivariate noncentral t distribution. This function follows Lan-DeMets for consistency with previous results.
}
\value{
  Drift value for the given condition
}
\examples{
  Drift(seqBound(ti=(1:5)/5)[, "up.bound"])
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
