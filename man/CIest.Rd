\name{CIest}
\alias{CIest}
\title{Confidence Interval Estimation}
\description{
  Get the point estimate and its confidence interval with a given contrast and alpha value using the t distribution.
}
\usage{
  CIest(Formula, Data, Term, Contrast, conf.level=0.95) 
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{Contrast}{a level vector. Levels are alphabetically ordered by default.}
  \item{conf.level}{confidence level of the confidence interval}
}
\details{
  Get the point estimate and its confidence interval with a given contrast and alpha value using the t distribution.
}
\value{
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for the t distribution}
  \item{Df}{degrees of freedom}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual degrees of freedom}
}
\examples{
  CIest(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, "TRT", c(-1, 1), 0.90) # 90% CI
}
\author{Kyun-Seop Bae k@acr.kr}
