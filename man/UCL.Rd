\name{UCL}
\alias{UCL}
\title{Upper Confidence Limit}
\description{
  The estimate of the upper bound of the confidence limit using t-distribution
}
\usage{
  UCL(y, conf.level=0.95)
}
\arguments{
  \item{y}{a vector of numerics}
  \item{conf.level}{confidence level}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  The estimate of the upper bound of the confidence limit using t-distribution
}
\author{Kyun-Seop Bae k@acr.kr}
