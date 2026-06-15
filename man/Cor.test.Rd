\name{Cor.test}
\alias{Cor.test}
\title{Correlation test of multiple numeric columns}
\description{
  Testing correlation between numeric columns of data with the Pearson method.
}
\usage{
  Cor.test(Data, conf.level=0.95)
}
\arguments{
  \item{Data}{a matrix or a \code{data.frame}}
  \item{conf.level}{confidence level}
}
\details{
  It uses all numeric columns of the input data. It uses "pairwise.complete.obs" rows.
}
\value{
  Row names show which columns are used for the test.
  \item{Estimate}{point estimate of correlation}
  \item{Lower CL}{lower confidence limit}
  \item{Upper CL}{upper confidence limit}
  \item{t value}{t value of the t distribution}
  \item{Df}{degrees of freedom}
  \item{Pr(>|t|)}{probability with the t distribution}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  Cor.test(mtcars)
}
