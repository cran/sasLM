\name{CONTR}
\alias{CONTR}
\title{F Test with a Set of Contrasts}
\description{
  Do F test with a given set of contrasts.
}
\usage{
  CONTR(L, Formula, Data, mu=0)
}
\arguments{
  \item{L}{contrast matrix. Each row is a contrast.}
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{mu}{a vector of mu for the hypothesis L. The length should be equal to the row count of L.}  
}
\details{
  It performs F test with a given set of contrasts (a matrix). It is similar to the CONTRAST clause of SAS PROC GLM. This can test the hypotheis that the linear combination (function)'s mean vector is mu.
}
\value{
  Returns sum of square and its F value and p-value.
  \item{Df}{degree of freedom}
  \item{Sum Sq}{sum of square for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{proability of larger than F value}
}
\examples{
  CONTR(t(c(0, -1, 1)), uptake ~ Type, CO2) # sum of square 
  ANOVA(uptake ~ Type, CO2) # compare with the above
}
\seealso{\code{\link{cSS}}}
\author{Kyun-Seop Bae k@acr.kr}
