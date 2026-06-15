\name{e2}
\alias{e2}
\title{Get a Contrast Matrix for Type II SS}
\description{
  Makes a contrast matrix for Type II SS.
}
\usage{
  e2(x, eps=1e-8)
}
\arguments{
  \item{x}{an output of \code{ModelMatrix}}
  \item{eps}{A value less than this is considered zero.}
}
\details{
  It makes a contrast matrix for Type II SS. If \code{zapsmall} is used, the result becomes less accurate.
}
\value{
  A contrast matrix for Type II SS.
}
\examples{
  round(e2(ModelMatrix(uptake ~ Plant + Type + Treatment + conc, CO2)), 12)
  round(e2(ModelMatrix(uptake ~ Type, CO2)), 12)
  round(e2(ModelMatrix(uptake ~ Type - 1, CO2)), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
