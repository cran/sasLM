\name{e2}
\alias{e2}
\title{Get a Contrast Matrix for Type II SS}
\description{
  Makes a contrast matrix for type II SS.
}
\usage{
  e2(x, eps=1e-8)
}
\arguments{
  \item{x}{an output of ModelMatrix}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It makes a contrast matrix for type II SS. If zapsmall is used, the result becomes more inaccurate.
}
\value{
  A contrast matrix for type II SS.
}
\examples{
  round(e2(ModelMatrix(uptake ~ Plant + Type + Treatment + conc, CO2)), 12)
  round(e2(ModelMatrix(uptake ~ Type, CO2)), 12)
  round(e2(ModelMatrix(uptake ~ Type - 1, CO2)), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
