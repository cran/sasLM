\name{REG}
\alias{REG}
\title{Regression of Linear Least Square, similar to SAS PROC REG}
\description{
  REG is similar to SAS PROC REG.
}
\usage{
REG(Formula, Data, NOINT=FALSE, eps=1e-8, summarize=TRUE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{NOINT}{If NOINT is TRUE, no intercept model will be used. Always -1 or +0 will be ignored in the formula.}
  \item{eps}{Less than this value is considered as zero.}
  \item{summarize}{If this is \code{FALSE}, \code{REG} returns just \code{lfit} result.}
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
  
  If \code{summarize=FALSE}, \code{REG} returns;
  \item{coeffcients}{beta coefficients}
  \item{g2}{g2 inverse}
  \item{rank}{rank of the model matrix}
  \item{DFr}{degree of freedom for the residual}
  \item{SSE}{sum of square error}
}
\examples{
  REG(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
