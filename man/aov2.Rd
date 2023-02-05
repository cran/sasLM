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
  \item{BETA}{if \code{TRUE}, coefficients (parameters) of \code{REG} will be returned. This is equivalent to SOLUTION option of SAS PROC GLM}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
}
\details{
  It performs the core function of SAS PROC ANOVA.
}
\value{
  The result table is comparable to that of SAS PROC ANOVA.
  \item{Df}{degree of freedom}
  \item{Sum Sq}{sum of square for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{proability of larger than F value} 

  Next returns are optional.
  \item{Parameter}{Parameter table with standard error, t value, p value. \code{TRUE} is 1, and \code{FALSE} is 0 in the \code{Estimable} column. This is returned only with BETA=TRUE option.}
  \item{Fitted}{Fitted value or y hat. This is returned only with Resid=TRUE option.}
  \item{Residual}{Weigthed residuals. This is returned only with Resid=TRUE option.}
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
