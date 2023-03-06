\name{trimmedMean}
\alias{trimmedMean}
\title{Trimmed Mean}
\description{
  Trimmed mean wrapping \code{mean} function.
}
\usage{
  trimmedMean(y, Trim=0.05)
}
\arguments{
  \item{y}{a vector of numerics}
  \item{Trim}{trimming proportion. Default is 0.05}
}
\details{
  It removes \code{NA} in the input vector.
}
\value{
  The value of trimmed mean
}
\author{Kyun-Seop Bae k@acr.kr}
