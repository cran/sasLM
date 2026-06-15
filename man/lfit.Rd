\name{lfit}
\alias{lfit}
\title{Linear Fit}
\description{
  Fits a least squares linear model.
}
\usage{
  lfit(x, y, eps=1e-8)
}
\arguments{
  \item{x}{a result of \code{ModelMatrix}}
  \item{y}{a column vector of the response (dependent) variable}
  \item{eps}{A value less than this is considered zero.}
}
\details{
  A minimal version of the least squares fit of a linear model
}
\value{
  \item{coefficients}{beta coefficients}
  \item{g2}{g2 inverse}
  \item{rank}{rank of the model matrix}
  \item{DFr}{degrees of freedom for the residual}
  \item{SSE}{sum of squares error}
  \item{SST}{sum of squares total}
  \item{DFr2}{degrees of freedom of the residual for the beta coefficient}
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{
  \code{\link{ModelMatrix}}
}
\examples{
  f1 = uptake ~ Type*Treatment + conc
  x = ModelMatrix(f1, CO2)
  y = model.frame(f1, CO2)[,1]
  lfit(x, y)
}
