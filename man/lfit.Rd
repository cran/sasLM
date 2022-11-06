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
  \item{SSE}{sum of squares error}
  \item{SST}{sum of squares total}
  \item{DFr2}{degree of freedom of the residual for beta coefficient}
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
