\name{regD}
\alias{regD}
\title{Regression of Conventional Way with Rich Diagnostics}
\description{
  \code{regD} provides rich diagnostics such as student residual, leverage(hat), Cook's D, studentized deleted residual, DFFITS, and DFBETAS.
}
\usage{
regD(Formula, Data)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
}
\details{
  It performs the conventional regression analysis. This does not use g2 inverse, therefore it cannot handle singular matrix. If the model(design) matrix is not full rank, use \code{REG} or less parameters.
}
\value{
  \item{Coefficients}{conventional coefficients summary with Wald statistics}
  \item{Diagnostics}{Diagnostics table for detecting outlier or influential/leverage points. This includes fitted (Predicted), residual (Residual), standard error of residual(se_resid), studentized residual(RStudent), hat(Leverage), Cook's D, studentized deleted residual(sdResid), DIFFITS, and COVRATIO.}
  \item{DFBETAS}{Column names are the names of coefficients. Each row shows how much each coefficient is affected by deleting the coressponding row of observation.}
}
\examples{
  regD(uptake ~ conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
