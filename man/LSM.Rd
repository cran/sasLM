\name{LSM}
\alias{LSM}
\title{Least Square Means}
\description{
  Estimates least square means using g2 inverse.
}
\usage{
  LSM(Formula, Data, Term, conf.level=0.95, adj="lsd", hideNonEst=TRUE, 
      PLOT=FALSE, descend=FALSE, ...)
}
\arguments{
  \item{Formula}{a conventional formula of model}
  \item{Data}{data.frame}
  \item{Term}{term name to be returned. If there is only one independent variable, this can be omitted.}
  \item{conf.level}{confidence level for the confidence limit}
  \item{adj}{adjustment method for grouping, "lsd"(default), "tukey", "bon", "duncan", "scheffe" are available. This does not affects SE, Lower CL, Upper CL of the output table.}  
  \item{hideNonEst}{logical. hide non-estimables}
  \item{PLOT}{logical. whether to plot LSMs and their confidence intervals}
	\item{descend}{logical. This specifies the plotting order be ascending or descending.}
  \item{...}{arguments to be passed to \code{plot}}
}
\details{
  It corresponds to SAS PROC GLM LSMEANS. The result of the second example below may be different from \code{emmeans}. This is because SAS or this function calculates mean of the transformed continuous variable. However, \code{emmeans} calculates the average before the transformation. Interaction of three variables is not supported yet. For adjust method "dunnett", see \code{PDIFF} function.
}
\value{
  Returns a table of expectations, t values and p-values.
  \item{Group}{group character. This appears with one-way ANOVA or \code{Term} or \code{adj} argument is provided.  }
  \item{LSmean}{point estimate of least square mean}
  \item{LowerCL}{lower confidence limit with the given confidence level by "lsd" method}
  \item{UpperCL}{upper confidence limit with the given confidence level by "lsd" method}
  \item{SE}{standard error of the point estimate}
  \item{Df}{degree of freedom of point estimate}
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
