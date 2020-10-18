\name{CIest}
\alias{CIest}
\title{Confidence Interval Estimation}
\description{
  Get point estimate and its confidence interval with given contrast and alpha value using t distribution.
}
\usage{
  CIest(Formula, Data, Term, Contrast=c(-1, 1), Alpha=0.10) 
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{Contrast}{a level vector. Level is alphabetically ordered by default.}
  \item{Alpha}{0.05 means 95 percent and 0.10 means 90 percent confidence interval.}
}
\details{
  Get point estimate and its confidence interval with given contrast and alpha value using t distribution.
}
\value{
  Returns point estimate and its confidence interval
}
\examples{
  CIest(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata, Term="TRT") # 90% CI
}
\author{Kyun-Seop Bae k@acr.kr}
