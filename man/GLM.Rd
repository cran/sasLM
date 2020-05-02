\name{GLM}
\alias{GLM}
\title{General Linear Model similar to SAS PROC GLM}
\description{
  GLM is the main function of this package.
}
\usage{
GLM(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It performs the core function of SAS PROC GLM.
}
\value{
  The result is comparable to that of SAS PROC GLM.
  \item{ANOVA}{ANOVA table for the model}
  \item{Type I}{Type I sum of square table}
  \item{Type II}{Type II sum of square table}
  \item{Type III}{Type III sum of square table}
  \item{Parameter}{Parameter table with standard error, t value, p value}
}
\examples{
  GLM(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
