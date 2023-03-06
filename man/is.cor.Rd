\name{is.cor}
\alias{is.cor}
\title{Is it a correlation matrix?}
\description{
  Testing if the input matrix is a correlation matrix or not
}
\usage{
  is.cor(m, eps=1e-16)
}
\arguments{
  \item{m}{a presumed correlation matrix}
  \item{eps}{epsilon value. An absolute value less than this is considered as zero.}
}
\details{
  A diagonal component should not be necessarily 1. But it should be close to 1.
}
\value{
  TRUE or FALSE
}
\author{Kyun-Seop Bae k@acr.kr}
