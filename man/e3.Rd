\name{e3}
\alias{e3}
\title{Get a Contrast Matrix for Type III SS}
\description{
  Makes a contrast matrix for Type III SS.
}
\usage{
  e3(x, eps=1e-8)
}
\arguments{
  \item{x}{an output of \code{ModelMatrix}}
  \item{eps}{A value less than this is considered zero.}
}
\details{
  It makes a contrast matrix for Type III SS. If \code{zapsmall} is used, the result becomes less accurate.
}
\value{
  A contrast matrix for Type III SS.
}
\examples{
  round(e3(ModelMatrix(uptake ~ Plant + Type + Treatment + conc, CO2)), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
