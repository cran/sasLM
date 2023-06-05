\name{ztest}
\alias{ztest}
\title{Test for the difference of two groups' means}
\description{
  This is similar to two groups t-test, but using standard normal (Z) distribution.
}
\usage{
  ztest(m1, s1, n1, m0, s0, n0, conf.level=0.95, nullHypo=0)
}
\arguments{
  \item{m1}{mean of the first (test, active, experimental) group}
  \item{s1}{known standard deviation of the first group}
  \item{n1}{sample size of the first group}
  \item{m0}{mean of the second (reference, control, placebo) group}
  \item{s0}{known standard deviationo of the second group}
  \item{n0}{sample size of the second group}
  \item{conf.level}{confidence level}
  \item{nullHypo}{value for the difference of means under null hypothesis}
}
\details{
  Use this only for known standard deviations (or variances) or very large sample sizes per group.
}
\value{
  The output format is very similar to t.test
}
\examples{
  ztest(5.4, 10.5, 3529, 5.1, 8.9, 5190) # NEJM 388;15 p1386
}
\seealso{\code{\link{mtest}}, \code{\link{tmtest}}, \code{\link{TTEST}}}
\keyword{t-test}
\author{Kyun-Seop Bae k@acr.kr}
