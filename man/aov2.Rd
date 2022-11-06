\name{aov2}
\alias{aov2}
\title{ANOVA with Type II SS}
\description{
  ANOVA with Type II SS.
}
\usage{
  aov2(Formula, Data)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
}
\details{
  It performs the core function of SAS PROC ANOVA.
}
\value{
  The result table is comparable to that of SAS PROC ANOVA.
  \item{Df}{degree of freedom}
  \item{Sum Sq}{sum of square for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{proability of larger than F value} 
}
\examples{
  aov2(uptake ~ Plant + Type + Treatment + conc, CO2)
  aov2(uptake ~ Type, CO2)
  aov2(uptake ~ Type - 1, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
