\name{e2}
\alias{e2}
\title{Get a Contrast Matrix for Type II SS}
\description{
  Makes a contrast matrix for type II SS.
}
\usage{
  e2(Formula, Data)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
}
\details{
  It makes a contrast matrix for type II SS.
}
\value{
  Returns a contrast matrix for type II SS.
}
\examples{
  round(e2(uptake ~ Plant + Type + Treatment + conc, CO2), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
