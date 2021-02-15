\name{UCL}
\alias{UCL}
\title{Upper Confidence Limit}
\description{
  The estimate of the upper bound of confidence limit using t-distribution
}
\usage{
  UCL(x, conf.level=0.95)
}
\arguments{
  \item{x}{a vector of numerics}
  \item{conf.level}{confidence level}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  The estimate of the upper bound of confidence limit using t-distribution
}
\author{Kyun-Seop Bae k@acr.kr}
