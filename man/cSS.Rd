\name{cSS}
\alias{cSS}
\title{Sum of Square with a Given Contrast Set}
\description{
  Calculates sum of squares of a contrast from a \code{lfit} result.
}
\usage{
  cSS(K, rx, mu=0, eps=1e-8)
}
\arguments{
  \item{K}{contrast matrix. Each row is a contrast.}
  \item{rx}{a result of \code{lfit} function}  
  \item{mu}{a vector of mu for the hypothesis K. The length should be equal to the row count of K.}  
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It calculates sum of squares with given a contrast matrix and a \code{lfit} result. It corresponds to SAS PROC GLM CONTRAST. This can test the hypotheis that the linear combination (function)'s mean vector is mu.
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
  rx = REG(uptake ~ Type, CO2, summarize=FALSE)
  cSS(t(c(0, -1, 1)), rx) # sum of square 
  ANOVA(uptake ~ Type, CO2) # compare with the above
}
\seealso{
  \code{\link{CONTR}}
}
\author{Kyun-Seop Bae k@acr.kr}
