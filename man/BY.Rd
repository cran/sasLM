\name{BY}
\alias{BY}
\title{Analysis BY variable}
\description{
  Functions such as GLM, REG, and aov1 can be run by levels of a variable.
}
\usage{
  BY(FUN, Formula, Data, By, ...)
}
\arguments{
  \item{FUN}{Function name to be called, such as GLM or REG}
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{By}{a variable name in the \code{Data}}
  \item{...}{arguments to be passed to \code{FUN} function}
}
\details{
  This mimics the BY clause of SAS procedures.
}
\value{
  a list of \code{FUN} function outputs. The names of the list are the levels of the \code{By} variable.
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  BY(GLM, uptake ~ Treatment + as.factor(conc), CO2, By="Type")
  BY(REG, uptake ~ conc, CO2, By="Type")
}
