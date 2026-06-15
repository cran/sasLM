\name{OBFBound}
\alias{OBFBound}
\title{O'Brien-Fleming bounds for cumulative Z-test in Group Sequential Design}
\description{
  Sequential O'Brien-Fleming upper bounds for the cumulative Z-test on accumulating data. Z values are correlated. This is usually used for group sequential design.
}
\usage{
  OBFBound(K, alpha=0.05, side=2, ti=NULL, c0=NULL)
}
\arguments{
  \item{K}{count of tests, including the final one}
  \item{alpha}{goal alpha value for the last test at time 0.}
  \item{side}{1=one-sided test, 2=two-sided test}
  \item{ti}{times for test. These should be in [0, 1]. If not specified, equal intervals are assumed.}
  \item{c0}{correlation matrix. If not specified, Brownian motion is assumed.}
}
\details{
  It calculates O'Brien-Fleming upper z-bounds and cumulative alpha-values for the repeated test in group sequential design.
}
\value{
  The result is a matrix.
  \item{ti}{time of test}
  \item{z}{O'Brien-Fleming upper z-bound}
  \item{cum.alpha}{cumulative alpha-value}
}
\examples{
  OBFBound(K=2)
  OBFBound(K=3)
  OBFBound(K=4)
  OBFBound(K=5)
}
\seealso{
  \code{\link{seqBound}}, \code{\link{PocockBound}}
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
