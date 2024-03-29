\name{GLM}
\alias{GLM}
\title{General Linear Model similar to SAS PROC GLM}
\description{
  GLM is the main function of this package.
}
\usage{
  GLM(Formula, Data, BETA=FALSE, EMEAN=FALSE, Resid=FALSE, conf.level=0.95,
      Weights=1)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{BETA}{if \code{TRUE}, coefficients (parameters) of \code{REG} will be returned. This is equivalent to SOLUTION option of SAS PROC GLM}
  \item{EMEAN}{if \code{TRUE}, least square means (or expected means) will be returned. This is equivalent to LSMEANS clause of SAS PROC GLM}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
  \item{conf.level}{confidence level for the confidence limit of the least square mean}
  \item{Weights}{weights for the weighted least square}
}
\details{
  It performs the core function of SAS PROC GLM. Least square means for the interaction term of three variables is not supported yet.
}
\value{
  The result is comparable to that of SAS PROC GLM.
  \item{ANOVA}{ANOVA table for the model}
  \item{Fitness}{Some measures of goodness of fit such as R-square and CV}
  \item{Type I}{Type I sum of square table}
  \item{Type II}{Type II sum of square table}
  \item{Type III}{Type III sum of square table}
  \item{Parameter}{Parameter table with standard error, t value, p value. \code{TRUE} is 1, and \code{FALSE} is 0 in the \code{Estimable} column. This is returned only with BETA=TRUE option.}
  \item{Expected Mean}{Least square (or expected) mean table with confidence limit. This is returned only with EMEAN=TRUE option.}
  \item{Fitted}{Fitted value or y hat. This is returned only with Resid=TRUE option.}
  \item{Residual}{Weigthed residuals. This is returned only with Resid=TRUE option.}
}
\examples{
GLM(uptake ~ Type*Treatment + conc, CO2[-1,]) # Making data unbalanced
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], EMEAN=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], Resid=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE, EMEAN=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE, EMEAN=TRUE, Resid=TRUE)
}
\author{Kyun-Seop Bae k@acr.kr}
