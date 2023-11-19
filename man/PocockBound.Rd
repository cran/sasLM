\name{PocockBound}
\alias{PocockBound}
\title{Pocock (fixed) Bound for the cumulative Z-test with a final target alpha-value}
\description{
  Cumulative alpha values with cumulative hypothesis test with a fixed upper bound z-value in group sequential design.
}
\usage{
  PocockBound(K=2, alpha=0.05, side=2)
}
\arguments{
  \item{K}{total number of tests}
  \item{alpha}{alpha value at the final test}
  \item{side}{1=one-side test, 2=two-side test}
}
\details{
  Pocock suggested a fixed upper bound z-value for the cumulative hypothesis test in group sequential designs.
}
\value{
  a fixed upper bound z-value for the K times repated hypothesis test with a final alpha-value.
  Attributes are;
  \item{ti}{time of test, Even-interval is assumed.}
  \item{cum.alpha}{cumulative alpha valued}
}
\examples{
  PocockBound(K=2) # Z-value of upper bound for the two-stage design
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
