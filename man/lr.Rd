\name{lr}
\alias{lr}
\title{Linear Regression with g2 inverse}
\description{
  Coefficients are calculated with the g2 inverse. The output is similar to \code{summary(lm())}.
}
\usage{
  lr(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{A value less than this is considered zero.}
}
\details{
  It uses \code{G2SWEEP} to get the g2 inverse. The result is similar to \code{summary(lm())} without options.
}
\value{
  The result is comparable to that of SAS PROC REG.
  \item{Estimate}{point estimate of parameters, coefficients}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for the t distribution}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual degrees of freedom}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  lr(uptake ~ Plant + Type + Treatment + conc, CO2)
  lr(uptake ~ Plant + Type + Treatment + conc - 1, CO2)
  lr(uptake ~ Type, CO2)
  lr(uptake ~ Type - 1, CO2)
}
