\name{mtest}
\alias{mtest}
\title{Independent two groups t-test similar to PROC TTEST with summarized input}
\description{
  This is comparable to SAS PROC TTEST except using summarized input (sufficient statistics).
}
\usage{
  mtest(m1, s1, n1, m0, s0, n0, conf.level=0.95)
}
\arguments{
  \item{m1}{mean of the first (test, active, experimental) group}
  \item{s1}{sample standard deviation of the first group}
  \item{n1}{sample size of the first group}
  \item{m0}{mean of the second (reference, control, placebo) group}
  \item{s0}{sample standard deviation of the second group}
  \item{n0}{sample size of the second group}
  \item{conf.level}{confidence level}
}
\details{
  This uses summarized input. This also produces confidence intervals of means and variances by group.
}
\value{
  The output format is comparable to SAS PROC TTEST.
}
\examples{
  mtest(5.4, 10.5, 3529, 5.1, 8.9, 5190) # NEJM 388;15 p1386
}
\seealso{\code{\link{TTEST}}, \code{\link{tmtest}}, \code{\link{ztest}}}
\keyword{t-test}
\author{Kyun-Seop Bae k@acr.kr}
