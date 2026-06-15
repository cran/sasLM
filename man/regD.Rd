\name{regD}
\alias{regD}
\title{Regression of Conventional Way with Rich Diagnostics}
\description{
  \code{regD} provides rich diagnostics such as student residual, leverage (hat), Cook's D, studentized deleted residual, DFFITS, and DFBETAS.
}
\usage{
regD(Formula, Data)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
}
\details{
  It performs the conventional regression analysis. This does not use the g2 inverse; therefore, it cannot handle a singular matrix. If the model (design) matrix is not of full rank, use \code{REG} or fewer parameters.
}
\value{
  \item{Coefficients}{conventional coefficients summary with Wald statistics}
  \item{Diagnostics}{Diagnostics table for detecting outliers or influential/leverage points. This includes the fitted value (Predicted), residual (Residual), standard error of the residual (SE_Resid), studentized residual (Student_Res), hat (Leverage), Cook's D, studentized deleted residual (RStudent), DFFITS, and COVRATIO.}
  \item{DFBETAS}{Column names are the names of coefficients. Each row shows how much each coefficient is affected by deleting the corresponding observation.}
}
\examples{
  regD(uptake ~ conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
