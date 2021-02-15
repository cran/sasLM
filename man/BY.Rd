\name{BY}
\alias{BY}
\title{Analysis BY variable}
\description{
  GLM, REG, aov1 etc. functions can be run by levels of a variable.
}
\usage{
BY(FUN, Formula, Data, By, ...)
}
\arguments{
  \item{FUN}{Function name to be called such as GLM, REG}
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{By}{a variable name in the \code{Data}}
  \item{...}{arguments to be passed to \code{FUN} function}
}
\details{
  This mimics SAS procedues' BY clause.
}
\value{
  a list of \code{FUN} function outputs. The names are after each level.
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
BY(GLM, uptake ~ Treatment + as.factor(conc), CO2, By="Type")
BY(REG, uptake ~ conc, CO2, By="Type")
}
