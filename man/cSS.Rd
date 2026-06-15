\name{cSS}
\alias{cSS}
\title{Sum of Square with a Given Contrast Set}
\description{
  Calculates the sum of squares of a contrast from an \code{lfit} result.
}
\usage{
  cSS(K, rx, mu=0, eps=1e-8)
}
\arguments{
  \item{K}{contrast matrix. Each row is a contrast.}
  \item{rx}{a result of the \code{lfit} function}  
  \item{mu}{a vector of mu for the hypothesis K. The length should be equal to the row count of K.}  
  \item{eps}{Values less than this are considered zero.}
}
\details{
  It calculates the sum of squares with a given contrast matrix and an \code{lfit} result. It corresponds to SAS PROC GLM CONTRAST. This can test the hypothesis that the linear combination (function)'s mean vector is mu.
}
\value{
  Returns the sum of squares and its F value and p-value.
  \item{Df}{degrees of freedom}
  \item{Sum Sq}{sum of squares for the set of contrasts}
  \item{Mean Sq}{mean square}
  \item{F value}{F value for the F distribution}
  \item{Pr(>F)}{probability of a larger F value}
}
\examples{
  rx = REG(uptake ~ Type, CO2, summarize=FALSE)
  cSS(t(c(0, -1, 1)), rx) # sum of square 
  GLM(uptake ~ Type, CO2) # compare with the above
}
\seealso{
  \code{\link{CONTR}}
}
\author{Kyun-Seop Bae k@acr.kr}
