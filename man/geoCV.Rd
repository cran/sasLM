\name{geoCV}
\alias{geoCV}
\title{Geometric Coefficient of Variation in percentage}
\description{
  Geometric coefficient of variation in percentage.
}
\usage{
  geoCV(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\details{
  It removes \code{NA}. This is sqrt(exp(var(log(x))) - 1)*100.
}
\value{
  Geometric coefficient of variation in percentage.  
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
CV(mtcars$mpg)
}
