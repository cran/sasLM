\name{geoCV}
\alias{geoCV}
\title{Geometric Coefficient of Variation in percentage}
\description{
  Geometric coefficient of variation in percentage.
}
\usage{
  geoCV(y)
}
\arguments{
  \item{y}{a numeric vector}
}
\details{
  It removes \code{NA}. This is sqrt(exp(var(log(x))) - 1)*100.
}
\value{
  Geometric coefficient of variation in percentage.  
}
\seealso{
  \code{\link{geoMean}}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  geoCV(mtcars$mpg)
}

