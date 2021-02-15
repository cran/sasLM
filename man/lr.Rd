\name{lr}
\alias{lr}
\title{Linear Regression with g2 inverse}
\description{
  Coefficients calculated with g2 inverse. Output is similar to \code{summary(lm())}.
}
\usage{
lr(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It uses G2SWEEP to get g2 inverse. The result is similar to \code{summary(lm())} without options.
}
\value{
  The result is comparable to that of SAS PROC REG.
  \item{Estimate}{point estimate of parameters, coefficients}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Pr(>|t|}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
lr(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
