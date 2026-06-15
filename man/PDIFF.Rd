\name{PDIFF}
\alias{PDIFF}
\title{Pairwise Difference}
\description{
  Estimates pairwise differences by a common method.
}
\usage{
  PDIFF(Formula, Data, Term, conf.level=0.95, adj="lsd", ref, PLOT=FALSE, 
        reverse=FALSE, ...)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{conf.level}{confidence level of the confidence interval}
  \item{adj}{"lsd", "tukey", "scheffe", "bon", "duncan", or "dunnett" to adjust the p-value and confidence limit}
  \item{ref}{reference or control level for the Dunnett test. If missing, the first level is used, as SAS does.}
  \item{PLOT}{whether to plot the diffogram}
  \item{reverse}{reverse A - B to B - A}
  \item{...}{arguments to be passed to \code{plot}}
}
\details{
  It corresponds to the PDIFF option of SAS PROC GLM.
}
\value{
  Returns a table of expectations, t values and p-values. Output columns may vary according to the adjustment option.
  \item{Estimate}{point estimate of the input linear contrast}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{Std. Error}{standard error of the point estimate}
  \item{t value}{value for the t distribution}
  \item{Df}{degrees of freedom}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual's degrees of freedom}
}
\examples{
  PDIFF(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)")
  PDIFF(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)", adj="tukey")
}
\seealso{
  \code{\link{LSM}}, \code{\link{Diffogram}}
}
\author{Kyun-Seop Bae k@acr.kr}
