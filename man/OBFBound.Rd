\name{OBFBound}
\alias{OBFBound}
\title{O'Brien-Flemming bounds for cumulative Z-test in Group Sequential Design}
\description{
  Sequential O'Brien-Flemming upper bounds for cumulative Z-test on accumaltive data. Z values are correlated. This is usually used for group sequential design.
}
\usage{
  OBFBound(K, alpha=0.05, side=2, ti=NULL, c0=NULL)
}
\arguments{
  \item{K}{count of tests including the final}
  \item{alpha}{goal alpha value for the last test at time 0.}
  \item{side}{1=one-side test, 2=two-side test}
  \item{ti}{times for test. These should be [0, 1]. If not specified, equal interval is assumed.}
  \item{c0}{correlation matrix. if not specified, Brownian motion assummed.}
}
\details{
  It calculates O'Brien-Flemming upper z-bounds and cumulative alpha-values for the repeated test in group sequential design.
}
\value{
  The result is a matrix.
  \item{ti}{time of test}
  \item{z}{O'Brien-Flemming upper z-bound}
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
