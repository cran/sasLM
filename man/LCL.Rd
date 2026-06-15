\name{LCL}
\alias{LCL}
\title{Lower Confidence Limit}
\description{
  The estimate of the lower bound of the confidence limit using the t-distribution
}
\usage{
  LCL(y, conf.level=0.95)
}
\arguments{
  \item{y}{a vector of numerics}
  \item{conf.level}{confidence level}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  The estimate of the lower bound of the confidence limit using the t-distribution
}
\seealso{
  \code{\link{UCL}}
}
\author{Kyun-Seop Bae k@acr.kr}
