\name{pB}
\alias{pB}
\title{Plot Confidence and Prediction Bands for Simple Linear Regression}
\description{
  It plots bands of the confidence interval and prediction interval for simple linear regression.
}
\usage{
  pB(Formula, Data, Resol=300, conf.level=0.95, lx, ly, ...)
}
\arguments{
  \item{Formula}{a formula}
  \item{Data}{a data.frame}
  \item{Resol}{resolution for the output}
  \item{conf.level}{confidence level}
  \item{lx}{x position of legend}
  \item{ly}{y position of legend}
  \item{...}{arguments to be passed to \code{plot}}
}
\details{
  It plots. Discard return values. If \code{lx} or \code{ly} is missing, the legend position is calculated automatically.
}
\value{
  Ignore return values.
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
pB(hp ~ disp, mtcars)
pB(mpg ~ disp, mtcars)
}