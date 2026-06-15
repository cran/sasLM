\name{LSM}
\alias{LSM}
\title{Least Squares Means}
\description{
  Estimates least squares means using the g2 inverse.
}
\usage{
  LSM(Formula, Data, Term, conf.level=0.95, adj="lsd", hideNonEst=TRUE, 
      PLOT=FALSE, descend=FALSE, ...)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a term name to be returned. If there is only one independent variable, this can be omitted.}
  \item{conf.level}{confidence level for the confidence limit}
  \item{adj}{adjustment method for grouping; "lsd" (default), "tukey", "bon", "duncan", and "scheffe" are available. This does not affect the SE, Lower CL, and Upper CL of the output table.}
  \item{hideNonEst}{logical. whether to hide non-estimable values}
  \item{PLOT}{logical. whether to plot LSMs and their confidence intervals}
	\item{descend}{logical. This specifies whether the plotting order is ascending or descending.}
  \item{...}{arguments to be passed to \code{plot}}
}
\details{
  It is equivalent to the LSMEANS statement of SAS PROC GLM. The result of the second example below may be different from \code{emmeans}. This is because SAS and this function calculate the mean of the transformed continuous variable. However, \code{emmeans} calculates the average before the transformation. An interaction of three variables is not supported yet. For the "dunnett" adjustment method, see the \code{\link{PDIFF}} function.
}
\value{
  Returns a table of expectations, t values and p-values.
  \item{Group}{group character. This appears when the model is one-way ANOVA or when the \code{Term} or \code{adj} argument is provided.}
  \item{LSmean}{point estimate of the least squares mean}
  \item{LowerCL}{lower confidence limit at the given confidence level by the "lsd" method}
  \item{UpperCL}{upper confidence limit at the given confidence level by the "lsd" method}
  \item{SE}{standard error of the point estimate}
  \item{Df}{degrees of freedom of the point estimate}
}
\seealso{
  \code{\link{PDIFF}}, \code{\link{Diffogram}}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  LSM(uptake ~ Type, CO2[-1,])
  LSM(uptake ~ Type - 1, CO2[-1,])
  LSM(uptake ~ Type*Treatment + conc, CO2[-1,])
  LSM(uptake ~ Type*Treatment + conc - 1, CO2[-1,])
  LSM(log(uptake) ~ Type*Treatment + log(conc), CO2[-1,])
  LSM(log(uptake) ~ Type*Treatment + log(conc) - 1, CO2[-1,])
  LSM(log(uptake) ~ Type*Treatment + as.factor(conc), CO2[-1,])
  LSM(log(uptake) ~ Type*Treatment + as.factor(conc) - 1, CO2[-1,])
  LSM(log(CMAX) ~ SEQ/SUBJ + PRD + TRT, BEdata)
  LSM(log(CMAX) ~ SEQ/SUBJ + PRD + TRT - 1, BEdata)
}
