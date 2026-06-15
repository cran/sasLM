\name{SS}
\alias{SS}
\title{Sum of Square}
\description{
  Sum of squares with ANOVA.
}
\usage{
  SS(x, rx, L, eps=1e-8)
}
\arguments{
  \item{x}{a result of \code{ModelMatrix} containing design information}
  \item{rx}{a result of \code{lfit}}
  \item{L}{linear hypothesis, a full matrix matching the information in \code{x}}
  \item{eps}{Values less than this are considered as zero.}
}
\details{
  It calculates the sum of squares and completes the ANOVA table.
}
\value{
  \item{ANOVA table}{a classical ANOVA table without the residual (Error) part.}
}
\seealso{
  \code{\link{ModelMatrix}},
  \code{\link{lfit}}
}
\author{Kyun-Seop Bae k@acr.kr}
