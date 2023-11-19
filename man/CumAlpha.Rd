\name{CumAlpha}
\alias{CumAlpha}
\title{Cumulative Alpha for the Fixed Z-value}
\description{
  Cumulative alpha values with repeated hypothesis with a fixed upper bound z-value.
}
\usage{
  CumAlpha(x, K=2, side=2)
}
\arguments{
  \item{x}{fixed upper z-value bound for the repeated hypothesis test}
  \item{K}{total number of tests}
  \item{side}{1=one-side test, 2=two-side test}
}
\details{
  It calculates cumulative alpha-values for the even-interval repeated hypothesis test with a fixed upper bound z-value. It assumes linear (proportional) increase of information amount and Brownian motion of z-value, i.e. the correlation is sqrt(t_i/t_j).
}
\value{
  The result is a matrix.
  \item{ti}{time of test, Even-interval is assumed.}
  \item{cum.alpha}{cumulative alpha valued}
}
\examples{
  CumAlpha(x=qnorm(1 - 0.05/2), K=10) # two-side Z-test with alpha=0.05 for ten times
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
