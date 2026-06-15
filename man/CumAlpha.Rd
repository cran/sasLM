\name{CumAlpha}
\alias{CumAlpha}
\title{Cumulative Alpha with various z's and ti}
\description{
  Cumulative alpha values for repeated hypothesis tests with changing bound z-values and times of test (ti).
}
\usage{
  CumAlpha(z, side=2, ti=NULL, c0=NULL, Seed=5)
}
\arguments{
  \item{z}{vector of upper z-value bounds for the repeated hypothesis test}
  \item{side}{1=one-sided test, 2=two-sided test}
  \item{ti}{vector of times (or information amount) of test. All values should be in [0, 1] and sorted. If not specified, equal intervals are assumed.}
  \item{c0}{correlation matrix. If not specified, Brownian motion is assumed.}
  \item{Seed}{seed value for the \code{mvtnorm::pmvnorm} function}
}
\details{
  It calculates cumulative alpha-values for the repeated hypothesis test with a vector of upper bound z-values. If the times of test are not specified, linear (proportional) increase of information amount and Brownian motion of z-values are assumed, i.e. the correlation is sqrt(t_i/t_j).
}
\value{
  The result is a matrix.
  \item{ti}{time of test}
  \item{cum.alpha}{cumulative alpha values}
}
\examples{
  CumAlpha(z=rep(qnorm(1 - 0.05/2), 10)) # two-side Z-test with alpha=0.05 for ten times
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
