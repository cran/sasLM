\name{Cor.test}
\alias{Cor.test}
\title{Correlation test of multiple numeric columns}
\description{
  Testing correlation between numeric columns of data with Pearson method.
}
\usage{
  Cor.test(Data, conf.level=0.95)
}
\arguments{
  \item{Data}{a matrix or a data.frame}
  \item{conf.level}{confidence level}
}
\details{
  It uses all numeric columns of input data. It uses "pairwise.complete.obs" rows.
}
\value{
  Row names show which columns are used for the test
  \item{Estimate}{point estimate of correlation}
  \item{Lower CL}{upper confidence limit}
  \item{Upper CL}{lower condidence limit}
  \item{t value}{t value of the t distribution}
  \item{Df}{degree of freedom}
  \item{Pr(>|t|)}{probability with the t distribution}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  Cor.test(mtcars)
}
