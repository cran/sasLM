\name{GLM}
\alias{GLM}
\title{General Linear Model similar to SAS PROC GLM}
\description{
  GLM is the main function of this package.
}
\usage{
GLM(Formula, Data, lsm=FALSE, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{lsm}{if \code{TRUE}, least square mean will be in the output}
  \item{conf.level}{confidence level for the confidence limit of the least square mean}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It performs the core function of SAS PROC GLM. Least square means for the tnteraction term of three variables is not supported yet.
}
\value{
  The result is comparable to that of SAS PROC GLM.
  \item{ANOVA}{ANOVA table for the model}
  \item{Type I}{Type I sum of square table}
  \item{Type II}{Type II sum of square table}
  \item{Type III}{Type III sum of square table}
  \item{Parameter}{Parameter table with standard error, t value, p value}
  \item{Least Square Mean}{Least square mean table with confindence limit}
}
\examples{
GLM(uptake ~ Type*Treatment + conc, CO2[-1,]) # Making data unbalanced
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], lsm=TRUE)
}
\author{Kyun-Seop Bae k@acr.kr}
