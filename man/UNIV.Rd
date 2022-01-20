\name{UNIV}
\alias{UNIV}
\title{Univariate Descriptive Statistics}
\description{
  Returns descriptive statistics of a numeric vector. 
}
\usage{
  UNIV(x, conf.level = 0.95)
}
\arguments{
  \item{x}{a numeric vector}
  \item{conf.level}{confidence level for confidence limit}
}
\details{
  A convenient and comprehensive descriptive statistics. NA is removed during the calculation. This is similar to SAS PROC UNIVARIATE.
}
\value{
  \item{nAll}{count of all element in the input vector}
  \item{nNA}{count of NA element}
  \item{nFinite}{count of finite numbers}
  \item{Mean}{mean excluding NA}
  \item{SD}{standard deviation excluding NA}
  \item{CV}{coefficient of variation in percent}
  \item{SEM}{standard error of the sample mean, sample mean divided by nFinite}
  \item{LowerConfLimit}{lower confidence limit of given confidence interval}
  \item{UpperConfLimit}{upper confidence limit of given confidence interval}
  \item{TrimmedMean}{trimmed mean with trimming 1 - confidence level}
  \item{Min}{minimum value}
  \item{Q1}{first quartile value}
  \item{Median}{median value}
  \item{Q3}{third quartile value}
  \item{Max}{maximum value}
  \item{Range}{range of finite numbers. maximum - minimum}
  \item{Skewness}{skewness}
  \item{SkewnessSE}{standard error of skewness}
  \item{Kurtosis}{kurtosis}
  \item{KurtosisSE}{kurtosis}
  \item{GeometricMean}{geometric mean, calculated only when all given values are positive.}
  \item{GeometricCV}{geometric coefficient of variation in percent, calculated only when all given values are positive.}
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  UNIV(lh)
}
\keyword{Descriptive Statistics}
