\name{aov3}
\alias{aov3}
\title{ANOVA with Type III SS}
\description{
  ANOVA with Type III SS.
}
\usage{
aov3(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
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
  aov3(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
