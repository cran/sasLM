\name{REG}
\alias{REG}
\title{Regression of Linear Least sqaure}
\description{
  REG is similar to SAS PROC REG.
}
\usage{
REG(Formula, Data, NOINT=FALSE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{NOINT}{If NOINT is TRUE, no intercept model will be used. Always -1 or +0 will be ignored in the formula.}
}
\details{
  It performs the core function of SAS PROC REG.
}
\value{
  The result is comparable to that of SAS PROC REG.
  \item{Estimate}{point estimate of parameters, coefficients}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Pr(>|t|}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
  REG(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
