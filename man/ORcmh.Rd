\name{ORcmh}
\alias{ORcmh}
\title{Odds Ratio of two groups with strata by CMH method}
\description{
  Odds ratio and its score confidence interval of two groups with stratification by Cochran-Mantel-Haenszel method
}
\usage{
  ORcmh(d0, conf.level=0.95)
}
\arguments{
  \item{d0}{A data.frame or matrix, of which each row means a strata. This should have four columns named y1, n1, y2, and n2; y1 and y2 for events of each group, n1 and n2 for sample size of each strata. The second group is usually the control group.}
  \item{conf.level}{confidence level}
}
\details{
  It calculates odds ratio and its score confidence interval of two groups. This can be used for meta-analysis also.
}
\value{
  The following output will be returned for each stratum and common value. There is no standard error.
  \item{odd1}{odd from the first group, y1/(n1 - y1)}
  \item{odd2}{odd from the second group, y2/(n2 - y2)}
  \item{OR}{odds ratio, odd1/odd2. The point estimate of common OR is calculated with MH weight.}
  \item{lower}{lower confidence limit of OR}
  \item{upper}{upper confidence limit of OR}
}
\examples{
  d1 = matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE)
  colnames(d1) =  c("y1", "n1", "y2", "n2")
  ORcmh(d1)
}
\seealso{
  \code{\link{RDmn1}}, \code{\link{RRmn1}}, \code{\link{ORmn1}}, \code{\link{RDmn}}, \code{\link{RRmn}}, \code{\link{ORmn}}, \code{\link{RDinv}}, \code{\link{RRinv}}, \code{\link{ORinv}}
}
\keyword{2x2 Table}
\author{Kyun-Seop Bae k@acr.kr}
