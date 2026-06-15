\name{ORmn}
\alias{ORmn}
\title{Odds Ratio and Score CI of two groups with strata by the MN method}
\description{
  Odds ratio and its score confidence interval of two groups with stratification by the Miettinen and Nurminen method
}
\usage{
  ORmn(d0, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{d0}{A \code{data.frame} or matrix, of which each row means a stratum. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each stratum. The second group is usually the control group.}
  \item{conf.level}{confidence level}
  \item{eps}{absolute value less than eps is regarded as negligible}
}
\details{
  It calculates the common odds ratio and its score confidence interval of two groups with stratification. The confidence interval is asymmetric, and there is no standard error in the output. For the stratified case, the inverse variance weighted score statistic with the bias correction is used, following Laud. The common odds ratio point estimate is the zero of the weighted score statistic, and the confidence bounds are found with the \code{uniroot} function. The result agrees with \code{ratesci::scoreci(contrast="OR", stratified=TRUE, skew=FALSE)} to at least 7 significant digits. For a single stratum, it returns the classical Miettinen-Nurminen interval of \code{\link{ORmn1}}. This can be used for meta-analysis also.
}
\value{
  The following output will be returned for each stratum and common value. There is no standard error.
  \item{odd1}{odds from the first group, y1/(n1 - y1). For the common value, it is calculated with the weights at the point estimate.}
  \item{odd2}{odds from the second group, y2/(n2 - y2). For the common value, it is calculated with the weights at the point estimate.}
  \item{OR}{odds ratio of the stratum. The common OR is the zero of the inverse variance weighted score statistic.}
  \item{lower}{lower confidence limit of OR}
  \item{upper}{upper confidence limit of OR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  ORmn(d1)
  d2 = data.frame(y1=c(4, 2, 10), n1=c(20, 20, 20), y2=c(8, 11, 2), n2=c(20, 20, 20))
  ORmn(d2)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{RDinv}}, \code{\link{RRinv}}, \code{\link{ORinv}}, \code{\link{ORcmh}}
}
\references{
  Miettinen O, Nurminen M. Comparative analysis of two rates. Stat Med 1985;4:213-26

  Laud PJ. Equal-tailed confidence intervals for comparison of rates. Pharmaceutical Statistics 2017;16:334-348
}
\keyword{2x2 Table}
\author{Kyun-Seop Bae k@acr.kr}
