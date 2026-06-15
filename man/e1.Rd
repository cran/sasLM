\name{e1}
\alias{e1}
\title{Get a Contrast Matrix for Type I SS}
\description{
  Makes a contrast matrix for Type I SS using the forward Doolittle method.
}
\usage{
  e1(XpX, eps=1e-8)
}
\arguments{
  \item{XpX}{the crossproduct of a design or model matrix. This should have appropriate column names.}
  \item{eps}{A value less than this is considered zero.}
}
\details{
  It makes a contrast matrix for Type I SS. If \code{zapsmall} is used, the result becomes less accurate.
}
\value{
  A contrast matrix for Type I SS.
}
\examples{
  x = ModelMatrix(uptake ~ Plant + Type + Treatment + conc, CO2)
  round(e1(crossprod(x$X)), 12)
}
\author{Kyun-Seop Bae k@acr.kr}
