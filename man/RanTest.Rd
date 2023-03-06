\name{RanTest}
\alias{RanTest}
\title{Test with Random Effects}
\description{
  Hypothesis test of with specified type SS using random effects as error terms. This corresponds to SAS PROC GLM's RANDOM /TEST clause.
}
\usage{
  RanTest(Formula, Data, Random="", Type=3, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Random}{a vector of random effects. All should be specified as primary terms, not as interaction terms. All interaction terms with random factor are regarded as random effects.}
  \item{Type}{Sum of square type to be used as contrast}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  Type can be from 1 to 3. All interaction terms with random factor are regarded as random effects.
  Here the error term should not be MSE.
}
\value{
  Returns ANOVA and E(MS) tables with specified type SS.
}
\examples{
  RanTest(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, Random="SUBJ")
  fBE = log(CMAX) ~ ADM/SEQ/SUBJ + PRD + TRT
  RanTest(fBE, BEdata, Random=c("ADM", "SUBJ"))
  RanTest(fBE, BEdata, Random=c("ADM", "SUBJ"), Type=2)
  RanTest(fBE, BEdata, Random=c("ADM", "SUBJ"), Type=1)
}
\author{Kyun-Seop Bae k@acr.kr}
