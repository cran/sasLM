\name{tmtest}
\alias{tmtest}
\title{Independent two means test similar to t.test with summarized input}
\description{
  This produces essentially the same to t.test except using summarized input (sufficient statistics).
}
\usage{
  tmtest(m1, s1, n1, m0, s0, n0, conf.level=0.95, nullHypo=0, var.equal=F)
}
\arguments{
  \item{m1}{mean of the first (test, active, experimental) group}
  \item{s1}{sample standard deviation of the first group}
  \item{n1}{sample size of the first group}
  \item{m0}{mean of the second (reference, control, placebo) group}
  \item{s0}{sample standard deviation of the second group}
  \item{n0}{sample size of the second group}
  \item{conf.level}{confidence level}
  \item{nullHypo}{value for the difference of means under null hypothesis}
  \item{var.equal}{assumption on the variance equality}
}
\details{
  The default is Welch t-test with Satterthwaite approximation.
}
\value{
  The output format is very similar to t.test
}
\examples{
  tmtest(5.4, 10.5, 3529, 5.1, 8.9, 5190) # NEJM 388;15 p1386
  tmtest(5.4, 10.5, 3529, 5.1, 8.9, 5190, var.equal=TRUE)
}
\seealso{\code{\link{mtest}}, \code{\link{TTEST}}, \code{\link{ztest}}}
\keyword{t-test}
\author{Kyun-Seop Bae k@acr.kr}
