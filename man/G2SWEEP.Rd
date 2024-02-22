\name{G2SWEEP}
\alias{G2SWEEP}
\title{Generalized inverse matrix of type 2 for linear regression}
\description{
  Generalized inverse is usually not unique. Some programs use this algorithm to get a unique generalized inverse matrix.
}
\usage{
  G2SWEEP(A, Augmented=FALSE, eps=1e-08) 
}
\arguments{
  \item{A}{a matrix to be inverted. If \code{A} is not a square matrix, G2SWEEP calls g2inv function.}
  \item{Augmented}{If this is \code{TRUE} and \code{A} is a model(design) matrix X, the last column should be X'y, the last row y'X, and the last cell y'y. See the reference and example for the detail. If the input matrix \code{A} is not a square matrix, \code{Augmented} option cannot be TRUE.}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  Generalized inverse of g2-type is used by some softwares to do linear regression. See 'SAS Technical Report R106, The Sweep Operator: Its importance in Statistical Computing' by J. H. Goodnight for the detail.
}
\value{
  \item{when Augmented=FALSE}{ordinary g2 inverse}
  \item{when Augmented=TRUE}{g2 inverse and beta hats in the last column and the last row, and sum of square error (SSE) in the last cell}
  \item{attribute "rank"}{the rank of input matrix}
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{
  \code{\link{lfit}}, \code{\link{ModelMatrix}}
}
\examples{
  f1 = uptake ~ Type + Treatment # formula
  x = ModelMatrix(f1, CO2)  # Model matrix and relevant information
  y = model.frame(f1, CO2)[, 1] # observation vector
  nc = ncol(x$X) # number of columns of model matrix
  XpY = crossprod(x$X, y)
  aXpX = rbind(cbind(crossprod(x$X), XpY), cbind(t(XpY), crossprod(y)))
  ag2 = G2SWEEP(aXpX, Augmented=TRUE)
  b = ag2[1:nc, (nc + 1)] ; b # Beta hat
  iXpX = ag2[1:nc, 1:nc] ; iXpX # g2 inverse of X'X
  SSE = ag2[(nc + 1), (nc + 1)] ; SSE # Sum of Square Error
  DFr = nrow(x$X) - attr(ag2, "rank") ; DFr # Degree of freedom for the residual

# Compare the below with the above
  REG(f1, CO2)
  aov1(f1, CO2)
}