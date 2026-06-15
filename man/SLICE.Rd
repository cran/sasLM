\name{SLICE}
\alias{SLICE}
\title{F Test with Slice}
\description{
  Performs an F test with a given slice term.
}
\usage{
  SLICE(Formula, Data, Term, By)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Term}{a factor name (not an interaction) to calculate the sum of squares and do an F test with least square means}
  \item{By}{a factor name to be used for slicing}
}
\details{
  It performs an F test with a given slice term. It is similar to the SLICE option of SAS PROC GLM.
}
\value{
  Returns the sum of squares and its F value and p-value. Row names are the levels of the slice term.
  \item{Df}{degrees of freedom}
  \item{Sum Sq}{sum of squares for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{probability of a larger F value}
}
\examples{
  SLICE(uptake ~ Type*Treatment, CO2, "Type", "Treatment") 
  SLICE(uptake ~ Type*Treatment, CO2, "Treatment", "Type") 
}
\author{Kyun-Seop Bae k@acr.kr}
