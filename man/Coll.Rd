\name{Coll}
\alias{Coll}
\title{Collinearity Diagnostics}
\description{
  Collearity digsnotics with tolerance, VIF, eigenvalue, condition index, variance proportions
}
\usage{
  Coll(Formula, Data)
}
\arguments{
  \item{Formula}{fomula of the model}
  \item{Data}{input data as a matrix or data.frame}
}
\details{
  Sometimes collinearity diagnostics after multiple linear regression are necessary.
}
\value{
  \item{Tol}{tolerance of independent variables}
  \item{VIF}{variance inflation factor of independent variables}
  \item{Eigenvalue}{eigenvalue of Z'Z (crossproduct) of standardized independent variables}
  \item{Cond. Index}{condition index}
  \item{Proportions of variances}{under the names of coefficients}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
Coll(mpg ~ disp + hp + drat + wt + qsec, mtcars) 
}
