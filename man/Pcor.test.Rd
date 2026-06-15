\name{Pcor.test}
\alias{Pcor.test}
\title{Partial Correlation test of multiple columns}
\description{
  Testing partial correlation between many columns of data with the Pearson method.
}
\usage{
  Pcor.test(Data, x, y)
}
\arguments{
  \item{Data}{a numeric matrix or \code{data.frame}}
  \item{x}{names of columns to be tested}
  \item{y}{names of control columns}
}
\details{
  It performs multiple partial correlation tests. It uses "complete.obs" rows of the x and y columns.
}
\value{
  Row names show which columns are used for the test.
  \item{Estimate}{point estimate of correlation}
  \item{Df}{degrees of freedom}
  \item{t value}{t value of the t distribution}
  \item{Pr(>|t|)}{probability with the t distribution}
}
\examples{
Pcor.test(mtcars, c("mpg", "hp", "qsec"), c("drat", "wt"))
}
\author{Kyun-Seop Bae k@acr.kr}
