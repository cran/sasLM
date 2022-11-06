\name{QuartileRange}
\alias{QuartileRange}
\title{Inter-Quartile Range}
\description{
  Interquartile range (Q3 - Q1) with a conventional formula.
}
\usage{
  QuartileRange(y, Type=2)
}
\arguments{
  \item{y}{a vector of numerics}
  \item{Type}{a type specifier to be passed to \code{IQR} function}
}
\details{
  It removes \code{NA} in the input vector. Type 2 is SAS default, while Type 6 is SPSS default.
}
\value{
  The value of interquartile range
}
\author{Kyun-Seop Bae k@acr.kr}
