\name{LSM}
\alias{LSM}
\title{Least Square Means}
\description{
  Estimates least square means using g2 inverse.
}
\usage{
  LSM(Formula, Data, conf.level=0.95, hideNonEst=TRUE)
}
\arguments{
  \item{Formula}{a conventional formula of model}
  \item{Data}{data.frame}
  \item{conf.level}{confidence level for the confidence limit}
  \item{hideNonEst}{hide non-estimables}
}
\details{
  It corresponds to SAS PROC GLM LSMEANS. The result of the second example below may be different from \code{emmeans}. This is because SAS or this function calculates mean of the transformed continuous variable. However, \code{emmeans} calculates the average before the transformation. Interaction of three variables is not supported yet.
}
\value{
  Returns a table of expectations, t values and p-values.
  \item{LSmean}{point estimate of least square mean}
  \item{LowerCL}{lower confidence limit with the given confidence level}
  \item{UpperCL}{upper confidence limit with the given confidence level}
  \item{SE}{standard error of the point estimate}
  \item{Df}{degree of freedom of point estimate}
}
\examples{
LSM(uptake ~ Type*Treatment + conc, CO2[-1,])
LSM(log(uptake) ~ Type*Treatment + log(conc), CO2[-1,])
LSM(log(uptake) ~ Type*Treatment + as.factor(conc), CO2[-1,])
LSM(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata)
}
\author{Kyun-Seop Bae k@acr.kr}
