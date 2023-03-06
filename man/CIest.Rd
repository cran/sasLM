\name{CIest}
\alias{CIest}
\title{Confidence Interval Estimation}
\description{
  Get point estimate and its confidence interval with given contrast and alpha value using t distribution.
}
\usage{
  CIest(Formula, Data, Term, Contrast, conf.level=0.95) 
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{Contrast}{a level vector. Level is alphabetically ordered by default.}
  \item{conf.level}{confidence level of confidence interval}
}
\details{
  Get point estimate and its confidence interval with given contrast and alpha value using t distribution.
}
\value{
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Df}{degree of freedom}
  \item{Pr(>|t|)}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
  CIest(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, "TRT", c(-1, 1), 0.90) # 90% CI
}
\author{Kyun-Seop Bae k@acr.kr}
