\name{T3test}
\alias{T3test}
\title{Test Type III SS using error term other than MSE}
\description{
  Hypothesis test of Type III SS using an error term other than MSE. This corresponds to SAS PROC GLM's RANDOM /TEST clause.
}
\usage{
  T3test(Formula, Data, Error="", eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Error}{an error term. Term name should be exactly same one listed the ANOVA output.}
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
  T3test(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, "SEQ:SUBJ")
}
\author{Kyun-Seop Bae k@acr.kr}
