\name{tsum1}
\alias{tsum1}
\title{Table Summary with 1 independent (x) variable}
\description{
  Summarize a continuous dependent (y) variable with one independent (x) variable.
}
\usage{
  tsum1(d, y, u, e=c("Mean", "SD", "N"), ou="", repl=list(c("length"), ("n")))
}
\arguments{
  \item{d}{a \code{data.frame} or matrix with column names}
  \item{y}{y variable name, a continuous variable}
  \item{u}{x variable name, upper side variable}
  \item{e}{a vector of summary function names}
  \item{ou}{order of levels of the upper side x variable}
  \item{repl}{a list of strings to replace after summarization. The length of the list should be 2, and both elements should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable with one x variable.
}
\value{
  A \code{data.frame} of summarized values. Row names are from \code{e} names. Column names are from the levels of the x variable.
}
\examples{
  tsum1(CO2, "uptake", "Treatment")
  tsum1(CO2, "uptake", "Treatment", 
        e=c("mean", "median", "sd", "min", "max", "length"), 
        ou=c("chilled", "nonchilled"),
        repl=list(c("median", "length"), c("med", "n")))
}
\seealso{\code{\link{tsum}}, \code{\link{tsum0}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\keyword{Table Summary}
\author{Kyun-Seop Bae k@acr.kr}
