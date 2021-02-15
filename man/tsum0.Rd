\name{tsum0}
\alias{tsum0}
\title{Table Summary 0 independent(x) variable}
\description{
  Summarize a continuous dependent(y) variable without any independent(x) variable.   
}
\usage{
  tsum0(d, y, e=c("mean", "sd", "length"), repl=list(c("length"), c("n")))
}
\arguments{
  \item{d}{a data.frame or matrix with colnames}
  \item{y}{y variable name, a continuous variable}
  \item{e}{a vector of summarize function names}
  \item{repl}{list of strings to replace after summarize. Length of list should be 2, and both should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable.
}
\value{
  A vector of summarized values 
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{\code{\link{tsum}}, \code{\link{tsum1}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\examples{
tsum0(CO2, "uptake")
tsum0(CO2, "uptake", repl=list(c("mean", "length"), c("Mean", "n")))
}
\keyword{Table Summary}

