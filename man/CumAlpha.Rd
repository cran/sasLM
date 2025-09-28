\name{CumAlpha}
\alias{CumAlpha}
\title{Cumulative Alpha with various z's and ti}
\description{
  Cumulative alpha values with repeated hypothesis with changing bound z-value and time of test (ti).
}
\usage{
  CumAlpha(z, side=2, ti=NULL, c0=NULL, Seed=5)
}
\arguments{
  \item{z}{vector of upper z-value bounds for the repeated hypothesis test}
  \item{side}{1=one-side test, 2=two-side test}
  \item{ti}{vector of times (or information amount) of test. All values should be [0, 1] and sorted. If not specified, equal interval is assumed.}
  \item{c0}{correlation matrix. if not specified, Brownian motion assummed.}
  \item{Seed}{seed value for mtvtnorm::pmvnorm function}
}
\details{
  It calculates cumulative alpha-values for the repeated hypothesis test with a vector of upper bound z-value. If the times of test is not specified, linear (proportional) increase of information amount and Brownian motion of z-value, i.e. the correlation is sqrt(t_i/t_j).
}
\value{
  The result is a matrix.
  \item{ti}{time of test}
  \item{cum.alpha}{cumulative alpha valued}
}
\examples{
  CumAlpha(z=rep(qnorm(1 - 0.05/2), 10)) # two-side Z-test with alpha=0.05 for ten times
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
