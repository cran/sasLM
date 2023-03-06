\name{Diffogram}
\alias{Diffogram}
\title{Plot Pairwise Differences}
\description{
  Plot pairwise differences by a common.
}
\usage{
  Diffogram(Formula, Data, Term, conf.level=0.95, adj="lsd", ...)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name to be estimated}
  \item{conf.level}{confidence level of confidence interval}
  \item{adj}{"lsd", "tukey", "scheffe", "bon", or "duncan" to adjust p-value and confidence limit}
  \item{...}{arguments to be passed to \code{plot}}
}
\details{
  This usually shows the shortest interval. It corresponds to SAS PROC GLM PDIFF. For adjust method "dunnett", see \code{PDIFF} function.
}
\value{
  no return value, but a plot on the current device
}
\examples{
  Diffogram(uptake ~ Type*Treatment + as.factor(conc), CO2, "as.factor(conc)")
}
\seealso{
  \code{\link{LSM}}, \code{\link{PDIFF}}
}
\author{Kyun-Seop Bae k@acr.kr}
