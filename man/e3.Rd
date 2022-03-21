\name{e3}
\alias{e3}
\title{Get a Contrast Matrix for Type III SS}
\description{
  Makes a contrast matrix for type III SS.
}
\usage{
  e3(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It makes a contrast matrix for type III SS. If zapsmall is used, the result becomes more inaccurate.
}
\value{
  A contrast matrix for type III SS.
}
\examples{
  round(e3(uptake ~ Plant + Type + Treatment + conc, CO2), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
