\name{SLICE}
\alias{SLICE}
\title{F Test with Slice}
\description{
  Do F test with a given slice term.
}
\usage{
  SLICE(Formula, Data, mTerm, sTerm)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{mTerm}{a factor name (not interaction) to calculate sum of square and do F test with least square means}
  \item{sTerm}{a factor name to be used for slice}  
}
\details{
  It performs F test with a given slice term. It is similar to the SLICE option SAS PROC GLM.
}
\value{
  Returns sum of square and its F value and p-value. Row names are the levels of the slice term.
  \item{Df}{degree of freedom}
  \item{Sum Sq}{sum of square for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{proability of larger than F value}
}
\examples{
  SLICE(uptake ~ Type*Treatment, CO2, "Type", "Treatment") 
  SLICE(uptake ~ Type*Treatment, CO2, "Treatment", "Type") 
}
\author{Kyun-Seop Bae k@acr.kr}
