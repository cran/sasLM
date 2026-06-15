\name{tsum0}
\alias{tsum0}
\title{Table Summary with 0 independent (x) variables}
\description{
  Summarize a continuous dependent (y) variable without any independent (x) variable.
}
\usage{
  tsum0(d, y, e=c("Mean", "SD", "N"), repl=list(c("length"), c("n")))
}
\arguments{
  \item{d}{a \code{data.frame} or matrix with column names}
  \item{y}{y variable name, a continuous variable}
  \item{e}{a vector of summary function names}
  \item{repl}{a list of strings to replace after summarization. The length of the list should be 2, and both elements should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable.
}
\value{
  A vector of summarized values 
}
\examples{
  tsum0(CO2, "uptake")
  tsum0(CO2, "uptake", repl=list(c("mean", "length"), c("Mean", "n")))
}
\seealso{\code{\link{tsum}}, \code{\link{tsum1}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\keyword{Table Summary}
\author{Kyun-Seop Bae k@acr.kr}
