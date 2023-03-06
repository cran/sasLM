\name{UNIV}
\alias{UNIV}
\title{Univariate Descriptive Statistics}
\description{
  Returns descriptive statistics of a numeric vector. 
}
\usage{
  UNIV(y, conf.level = 0.95)
}
\arguments{
  \item{y}{a numeric vector}
  \item{conf.level}{confidence level for confidence limit}
}
\details{
  A convenient and comprehensive function for descriptive statistics. NA is removed during the calculation. This is similar to SAS PROC UNIVARIATE.
}
\value{
  \item{nAll}{count of all elements in the input vector}
  \item{nNA}{count of NA element}
  \item{nFinite}{count of finite numbers}
  \item{Mean}{mean excluding NA}
  \item{SD}{standard deviation excluding NA}
  \item{CV}{coefficient of variation in percent}
  \item{SEM}{standard error of the sample mean, the sample mean divided by nFinite}
  \item{LowerConfLimit}{lower confidence limit of given confidence interval}
  \item{UpperConfLimit}{upper confidence limit of given confidence interval}
  \item{TrimmedMean}{trimmed mean with trimming 1 - confidence level}
  \item{Min}{minimum value}
  \item{Q1}{first quartile value}
  \item{Median}{median value}
  \item{Q3}{third quartile value}
  \item{Max}{maximum value}
  \item{Range}{range of finite numbers. maximum - minimum}
  \item{IQR}{inter-quartile range type 2, which is SAS default}
  \item{MAD}{mean absolute deviation}
  \item{Skewness}{skewness}
  \item{SkewnessSE}{standard error of skewness}
  \item{Kurtosis}{kurtosis}
  \item{KurtosisSE}{kurtosis}
  \item{GeometricMean}{geometric mean, calculated only when all given values are positive.}
  \item{GeometricCV}{geometric coefficient of variation in percent, calculated only when all given values are positive.}
}
\author{Kyun-Seop Bae k@acr.kr}
\keyword{Descriptive Statistics}
\examples{
  UNIV(lh)
}
