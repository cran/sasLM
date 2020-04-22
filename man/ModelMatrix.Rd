\name{ModelMatrix}
\alias{ModelMatrix}
\title{Model Matrix}
\description{
  This model matrix is similar to \code{model.matrix}. But it does not omit unnecessary columns.
}
\usage{
  ModelMatrix(Formula, Data, NOINT=FALSE, KeepOrder=FALSE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{NOINT}{If NOINT is TRUE, no intercept model will be used. Always -1 or +0 will be ignored in the formula.}
  \item{KeepOrder}{If \code{KeepOrder} is \code{TRUE}, terms in \code{Formula} will be kept. This is for Type I SS.}
}
\details{
  It makes the model(design) matrix for \code{GLM}.
}
\value{
  Model matrix and attributes similar to the output of \code{model.matrix}.
  \item{X}{design matrix, i.e. model matrix}
  \item{terms}{detailed information about terms such as formula and labels}
  \item{termsIndices}{term indices}
  \item{assign}{assignemnt of columns for each terms in order, different way of expressing term indices}
}
\author{Kyun-Seop Bae k@acr.kr}
