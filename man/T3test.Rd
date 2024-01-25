\name{T3test}
\alias{T3test}
\title{Test Type III SS using error term other than MSE}
\description{
  Hypothesis test of Type III SS using an error term other than MSE. This corresponds to SAS PROC GLM's RANDOM /TEST clause.
}
\usage{
  T3test(Formula, Data, H="", E="", eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{H}{Hypothesis term}
  \item{E}{Error term}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It tests a factor of type III SS using some other term as an error term.
  Here the error term should not be MSE.
}
\value{
  Returns one or more ANOVA table(s) of type III SS.
}
\examples{
  T3test(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, E=c("SEQ:SUBJ"))
  T3test(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, H="SEQ", E=c("SEQ:SUBJ"))
}
\author{Kyun-Seop Bae k@acr.kr}
