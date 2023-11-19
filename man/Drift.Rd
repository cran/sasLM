\name{Drift}
\alias{Drift}
\title{Drift defined by Lan and DeMets for Group Sequential Design}
\description{
  Calculate the drift value with given upper bounds (z-valuse), times of test, and power.
}
\usage{
  Drift(bi, ti=NULL, Power=0.9)
}
\arguments{
  \item{bi}{upper bound z-values}
  \item{ti}{times of test. These should be in the range of [0, 1]. If omitted, even-interval is assumed.}
  \item{Power}{target power at the final test}
}
\details{
  It calculates the drift value with given upper bound z-values, times of test, and power. If the times of test is not given, even-interval is assumed. \code{mvtnorm::pmvt} (with noncentrality) is better than pmvnorm in calculating power and sample size. But, Lan-DeMets used multi-variate normal rather than multi-variate noncentral t distributionh. This function followed Lan-DeMets for the consistency with previous results.
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
