\name{Pcor.test}
\alias{Pcor.test}
\title{Partial Correlation test of multiple columns}
\description{
  Testing partial correlation between many columns of data with Pearson method.
}
\usage{
  Pcor.test(Data, x, y)
}
\arguments{
  \item{Data}{a numeric matrix or data.frame}
  \item{x}{names of to be tested columns}
  \item{y}{names of control columns}
}
\details{
  It performs multiple partial correlation test. It uses "complete.obs" rows of x and y columns.
}
\value{
  Row names show which columns are used for the test
  \item{Estimate}{point estimate of correlation}
  \item{Df}{degree of freedom}
  \item{t value}{t value of the t distribution}
  \item{Pr(>|t|)}{probability with the t distribution}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
Pcor.test(mtcars, c("mpg", "hp", "qsec"), c("drat", "wt"))
}
