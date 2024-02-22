\name{g2inv}
\alias{g2inv}
\title{Generalized type 2 inverse matrix, g2 inverse}
\description{
  Generalized inverse is usually not unique. Some programs use this algorithm to get a unique generalized inverse matrix. This uses SWEEP operator and works for non-square matrix also.
}
\usage{
  g2inv(A, eps=1e-08) 
}
\arguments{
  \item{A}{a matrix to be inverted}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  See 'SAS Technical Report R106, The Sweep Operator: Its importance in Statistical Computing' by J. H. Goodnight for the detail.
}
\value{
  g2 inverse
}
\references{
  Searle SR, Khuri AI. Matrix Algebra Useful for Statistics. 2e. John Wiley and Sons Inc. 2017.
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{
  \code{\link{G2SWEEP}}
}
\examples{
  A = matrix(c(1, 2, 4, 3, 3, -1, 2, -2, 5, -4, 0, -7), byrow=TRUE, ncol=4) ; A
  g2inv(A)
}
