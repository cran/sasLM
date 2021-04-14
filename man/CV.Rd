\name{CV}
\alias{CV}
\title{Coefficient of Variation in percentage}
\description{
  Coefficient of variation in percentage.
}
\usage{
  CV(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\details{
  It removes \code{NA}.
}
\value{
  Coefficient of variation in percentage.  
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
CV(mtcars$mpg)
}
