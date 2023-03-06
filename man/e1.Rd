\name{e1}
\alias{e1}
\title{Get a Contrast Matrix for Type I SS}
\description{
  Makes a contrast matrix for type I SS using forward Doolittle method.
}
\usage{
  e1(XpX, eps=1e-8)
}
\arguments{
  \item{XpX}{crossprodut of a design or model matrix. This should have appropriate column names.}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It makes a contrast matrix for type I SS. If zapsmall is used, the result becomes more inaccurate.
}
\value{
  A contrast matrix for type I SS.
}
\examples{
  x = ModelMatrix(uptake ~ Plant + Type + Treatment + conc, CO2)
  round(e1(crossprod(x$X)), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
