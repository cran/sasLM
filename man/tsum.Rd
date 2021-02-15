\name{tsum}
\alias{tsum}
\title{Table Summary}
\description{
  Summarize a continuous dependent variable with or without independent variables.   
}
\usage{
  tsum(Formula=NULL, Data=NULL, ColNames=NULL, MaxLevel=30, ...)
}
\arguments{
  \item{Formula}{a conventional formula}
  \item{Data}{a data.frame or a matrix}
  \item{ColNames}{If there is no Formula, this will be used.}
  \item{MaxLevel}{More than this will not be handled.}
  \item{...}{arguments to be passed to \code{tsum0}, \code{tsum1}, \code{tsum2}, or \code{tsum3}}
}
\details{
  A convenient summarization function for a continuous variable. This is a wrapper function to \code{tsum0}, \code{tsum1}, \code{tsum2}, or \code{tsum3}.
}
\value{
  A data.frame of descriptive summarization values. 
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{\code{\link{tsum0}}, \code{\link{tsum1}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\examples{
tsum(lh)
t(tsum(CO2))
t(tsum(uptake ~ Treatment, CO2))
tsum(uptake ~ Type + Treatment, CO2)
print(tsum(uptake ~ conc + Type + Treatment, CO2), digits=3)
}
\keyword{Table Summary}
