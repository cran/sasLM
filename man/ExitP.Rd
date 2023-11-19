\name{ExitP}
\alias{ExitP}
\title{Exit Probability with cumulative Z-test in Group Sequential Design}
\description{
  Exit probabilities with given drift, upper bounds, and times of test.
}
\usage{
  ExitP(Theta, bi, ti=NULL)
}
\arguments{
  \item{Theta}{drift value defined by Lan-DeMets. See the reference.}
  \item{bi}{upper bound z-values}
  \item{ti}{times of test. These should be in the range of [0, 1]. If omitted, even-interval is assumed.}
}
\details{
  It calculates exit proabilities and cumulative exit probabilities with given drift, upper z-bounds and times of test. If the times of test is not given, even-interval is assumed. \code{mvtnorm::pmvt} (with noncentrality) is better than pmvnorm in calculating power and sample size. But, Lan-DeMets used multi-variate normal rather than multi-variate noncentral t distributionh. This function followed Lan-DeMets for the consistency with previous results.
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
