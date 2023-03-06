\name{ModelMatrix}
\alias{ModelMatrix}
\title{Model Matrix}
\description{
  This model matrix is similar to \code{model.matrix}. But it does not omit unnecessary columns.
}
\usage{
  ModelMatrix(Formula, Data, KeepOrder=FALSE, XpX=FALSE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{KeepOrder}{If \code{KeepOrder} is \code{TRUE}, terms in \code{Formula} will be kept. This is for Type I SS.}
  \item{XpX}{If \code{XpX} is \code{TRUE}, the cross-product of the design matrix (XpX, X'X) will be returned instead of the design matrix (X).}
}
\details{
  It makes the model(design) matrix for \code{GLM}.
}
\value{
  Model matrix and attributes similar to the output of \code{model.matrix}.
  \item{X}{design matrix, i.e. model matrix}
  \item{XpX}{cross-product of the design matrix, X'X}
  \item{terms}{detailed information about terms such as formula and labels}
  \item{termsIndices}{term indices}
  \item{assign}{assignemnt of columns for each term in order, different way of expressing term indices}
}
\author{Kyun-Seop Bae k@acr.kr}
