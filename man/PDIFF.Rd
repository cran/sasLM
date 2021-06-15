\name{PDIFF}
\alias{PDIFF}
\title{Pairwise Difference}
\description{
  Estimates pairwise difference by a common method.
}
\usage{
  PDIFF(Formula, Data, Term, conf.level=0.95, adj="lsd", ref, PLOT=FALSE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{conf.level}{confidence level of confidence interval}
  \item{adj}{"lsd", "tukey", "scheffe", "bon", "duncan", or "dunnett" to adjust p-value and confidence limit}  
  \item{ref}{reference or control level for Dunnett test}  
  \item{PLOT}{whether to plot or not the diffogram}  
}
\details{
  It corresponds to PDIFF option of SAS PROC GLM.
}
\value{
  Returns a table of expectations, t values and p-values. Outpuc columns may vary according to the adjustment option.
  \item{Estimate}{point estimate of the input linear constrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for t distribution}
  \item{Df}{degree of freedom}
  \item{Pr(>|t|)}{probability of larger than absolute t value from t distribution with residual's degree of freedom}
}
\examples{
  PDIFF(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)")
  PDIFF(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)", adj="tukey")
}
\author{Kyun-Seop Bae k@acr.kr}
