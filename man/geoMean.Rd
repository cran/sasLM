\name{geoMean}
\alias{geoMean}
\title{Geometric Mean without NA}
\description{
  mean without \code{NA} values.
}
\usage{
  geoMean(y)
}
\arguments{
  \item{y}{a vector of numerics}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  geometric mean value
}
\seealso{
  \code{\link{geoCV}}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  geoMean(mtcars$mpg)
}