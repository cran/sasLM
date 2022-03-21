\name{e2}
\alias{e2}
\title{Get a Contrast Matrix for Type II SS}
\description{
  Makes a contrast matrix for type II SS.
}
\usage{
  e2(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It makes a contrast matrix for type II SS. If zapsmall is used, the result becomes more inaccurate.
}
\value{
  A contrast matrix for type II SS.
}
\examples{
  round(e2(uptake ~ Plant + Type + Treatment + conc, CO2), 12)
  round(e2(uptake ~ Type, CO2), 12)
  round(e2(uptake ~ Type - 1, CO2), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
