\name{QuartileRange}
\alias{QuartileRange}
\title{Inter-Quartile Range}
\description{
  Interquartile range (Q3 - Q1) with a conventional formula.
}
\usage{
  QuartileRange(x, Type=6)
}
\arguments{
  \item{x}{a vector of numerics}
  \item{Type}{a type specifier to be passed to \code{IQR} function}
}
\details{
  It removes \code{NA} in the input vector. 
}
\value{
  The value of interquartile range
}
\author{Kyun-Seop Bae k@acr.kr}
