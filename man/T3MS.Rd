\name{T3MS}
\alias{T3MS}
\title{Type III Expected Mean Square Formula}
\description{
  Calculates a formula table for expected mean square of Type III SS.
}
\usage{
  T3MS(Formula, Data, L0, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{L0}{a matrix of row linear contrasts, if missed, \code{e3} is used}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  This is necessary for further hypothesis test of nesting factors.
}
\value{
  A coefficient matrix for Type III expected mean square
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  T3MS(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata)
}
