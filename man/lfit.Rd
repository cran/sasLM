\name{lfit}
\alias{lfit}
\title{Linear Fit}
\description{
  Fits a least square linear model. 
}
\usage{
  lfit(x, y, eps=1e-8)
}
\arguments{
  \item{x}{a result of ModelMatrix}
  \item{y}{a column vector of response, dependent variable}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  Minimum version of least square fit of a linear model
}
\value{
  \item{coeffcients}{beta coefficients}
  \item{g2}{g2 inverse}
  \item{rank}{rank of the model matrix}
  \item{DFr}{degree of freedom for the residual}
  \item{SSE}{sum of square error}
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{
  \code{\link{ModelMatrix}}
}