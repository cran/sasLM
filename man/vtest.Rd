\name{vtest}
\alias{vtest}
\title{F-Test for the ratio of two groups' variances}
\description{
  F-test for the ratio of two groups' variances. This is similar to var.test except using the summarized input.
}
\usage{
  vtest(v1, n1, v0, n0, ratio=1, conf.level=0.95)
}
\arguments{
  \item{v1}{sample variance of the first (test, active, experimental) group}
  \item{n1}{sample size of the first group}
  \item{v0}{sample variance of the second (reference, control, placebo) group}
  \item{n0}{sample size of the second group}
  \item{ratio}{value for the ratio of variances under null hypothesis}
  \item{conf.level}{confidence level}
}
\details{
  For the confidence interval of one group, use UNIV function.
}
\value{
  The output format is very similar to var.test.
}
\examples{
  vtest(10.5^2, 5190, 8.9^2, 3529) # NEJM 388;15 p1386
  vtest(2.3^2, 13, 1.5^2, 11, conf.level=0.9) # Red book p240
}
\author{Kyun-Seop Bae k@acr.kr}
