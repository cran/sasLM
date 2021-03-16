\name{estmb}
\alias{estmb}
\title{Estimability Check}
\description{
  Check the estimability of row vectors of coefficients.
}
\usage{
  estmb(L, X, g2, eps=1e-8)
}
\arguments{
  \item{L}{row vectors of coefficients}
  \item{X}{a model (design) matrix from \code{ModelMatrix}}
  \item{g2}{g2 generalized inverse of \code{crossprod(X)}}
  \item{eps}{absolute value less than this is considered to be zero.}
}
\details{
  It checks estimability of L, row vectors of coefficients. This corresponds to SAS PROC GLM ESTIMATE.
}
\value{
  a vector of logical values indicating which row is estimable (as TRUE)
}
\seealso{\code{\link{G2SWEEP}}}
\author{Kyun-Seop Bae k@acr.kr}
