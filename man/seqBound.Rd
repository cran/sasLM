\name{seqBound}
\alias{seqBound}
\title{Sequential bounds for cumulative Z-test in Group Sequential Design}
\description{
  Sequential upper bounds for cumulative Z-test on accumulating data. Z values are correlated. This is usually used for a group sequential design.
}
\usage{
  seqBound(ti, alpha = 0.05, side = 2, t2 = NULL, asf = 1)
}
\arguments{
  \item{ti}{times for the tests. These should be in [0, 1].}
  \item{alpha}{goal alpha value for the last test at time 1.}
  \item{side}{1=one-sided test, 2=two-sided test}
  \item{t2}{fractions of the information amount. These should be in [0, 1]. If not available, ti will be used instead.}
  \item{asf}{alpha spending function. 1=O'Brien-Fleming type (approximate, not exact), 2=Pocock type (approximate, not exact), 3=alpha*ti, 4=alpha*ti^1.5, 5=alpha*ti^2}
}
\details{
  It calculates upper z-bounds and cumulative alpha-values for the repeated tests in a group sequential design. The correlation is assumed to be sqrt(t_i/t_j).
  Use \code{\link{PocockBound}} and \code{\link{OBFBound}} for more exact bounds.
}
\value{
  The result is a matrix.
  \item{time}{time of test}
  \item{up.bound}{upper z-bound}
  \item{cum.alpha}{cumulative alpha-value}
}
\examples{
  seqBound(ti=(1:5)/5)
  seqBound(ti=(1:5)/5, asf=2)
}
\seealso{
  \code{\link{PocockBound}}, \code{\link{OBFBound}}
}
\references{
  Reboussin DM, DeMets DL, Kim K, Lan KKG. Computations for group sequential boundaries using the Lan-DeMets function method. Controlled Clinical Trials. 2000;21:190-207.
}
\keyword{Group Sequential Design}
\author{Kyun-Seop Bae k@acr.kr}
