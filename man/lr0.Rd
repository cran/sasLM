\name{lr0}
\alias{lr0}
\title{Simple Linear Regressions with Each Independent Variable}
\description{
  Usually, the first step to multiple linear regression is simple linear regressions with a single independent variable.
}
\usage{
  lr0(Formula, Data)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model. Intercept will always be added.}
  \item{Data}{a \code{data.frame} to be analyzed}
}
\details{
  It performs simple linear regression for each independent variable.
}
\value{
  Each row means one simple linear regression with that row name as the only independent variable.
  \item{Intercept}{estimate of the intecept}
  \item{SE(Intercept)}{standard error of the intercept}
  \item{Slope}{estimate of the slope}
  \item{SE(Slope)}{standard error of the slope}
  \item{Rsq}{R-squared for the simple linear model}
  \item{Pr(>F)}{p-value of slope or the model}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  lr0(uptake ~ Plant + Type + Treatment + conc, CO2)
  lr0(mpg ~ ., mtcars)
}
