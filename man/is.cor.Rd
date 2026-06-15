\name{is.cor}
\alias{is.cor}
\title{Is it a correlation matrix?}
\description{
  Tests if the input matrix is a correlation matrix or not.
}
\usage{
  is.cor(m, eps=1e-16)
}
\arguments{
  \item{m}{a presumed correlation matrix}
  \item{eps}{epsilon value. An absolute value less than this is considered zero.}
}
\details{
  A diagonal component does not need to be exactly 1, but it should be close to 1.
}
\value{
  TRUE or FALSE
}
\author{Kyun-Seop Bae k@acr.kr}
