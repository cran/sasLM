\name{aov2}
\alias{aov2}
\title{ANOVA with Type II SS}
\description{
  ANOVA with Type II SS.
}
\usage{
  aov2(Formula, Data, BETA=FALSE, Resid=FALSE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{BETA}{if \code{TRUE}, coefficients (parameters) of \code{REG} will be returned. This is equivalent to the SOLUTION option of SAS PROC GLM.}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
}
\details{
  It performs the core function of SAS PROC GLM, and returns Type II SS. This also accepts continuous independent variables.
}
\value{
  The result table is comparable to that of SAS PROC ANOVA.
  \item{Df}{degrees of freedom}
  \item{Sum Sq}{sum of squares for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{probability of a larger F value} 

  The next returns are optional.
  \item{Parameter}{Parameter table with standard error, t value, p value. \code{TRUE} is 1, and \code{FALSE} is 0 in the \code{Estimable} column. This is returned only with the \code{BETA=TRUE} option.}
  \item{Fitted}{Fitted values or y hat. This is returned only with the \code{Resid=TRUE} option.}
  \item{Residual}{Weighted residuals. This is returned only with the \code{Resid=TRUE} option.}
}
\examples{
  aov2(uptake ~ Plant + Type + Treatment + conc, CO2)
  aov2(uptake ~ Plant + Type + Treatment + conc, CO2, BETA=TRUE)
  aov2(uptake ~ Plant + Type + Treatment + conc, CO2, Resid=TRUE)
  aov2(uptake ~ Plant + Type + Treatment + conc, CO2, BETA=TRUE, Resid=TRUE)
  aov2(uptake ~ Type, CO2)
  aov2(uptake ~ Type - 1, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
