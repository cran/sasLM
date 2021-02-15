\name{PDIFF}
\alias{PDIFF}
\title{Pairwise Difference by Least Significant Difference}
\description{
  Estimates pairwise difference by least signicicant difference.
}
\usage{
  PDIFF(Formula, Data, Term, conf.level=0.95)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{conf.level}{confidence level of confidence interval}
}
\details{
  This usually shows the shortest interval. It corresponds to SAS PROC GLM PDIFF. Interaction of three variables is not supported yet.
}
\value{
  Returns a table of expectations, t values and p-values.
  \item{Estimate}{point estimate of the input linear constrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Df}{degree of freedom}
  \item{Pr(>|t|}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
PDIFF(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)")
}
\author{Kyun-Seop Bae k@acr.kr}
