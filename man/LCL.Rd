\name{LCL}
\alias{LCL}
\title{Lower Confidence Limit}
\description{
  The estimate of the lower bound of confidence limit using t-distribution
}
\usage{
  LCL(x, conf.level=0.95)
}
\arguments{
  \item{x}{a vector of numerics}
  \item{conf.level}{confidence level}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  The estimate of the lower bound of confidence limit using t-distribution
}
\author{Kyun-Seop Bae k@acr.kr}
