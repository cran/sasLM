\name{tsum1}
\alias{tsum1}
\title{Table Summary 1 independent(x) variable}
\description{
  Summarize a continuous dependent(y) variable with one independent(x) variable.   
}
\usage{
  tsum1(d, y, u, e=c("mean", "sd", "length"), ou="", 
        repl=list(c("length"), ("n")))
}
\arguments{
  \item{d}{a data.frame or matrix with colnames}
  \item{y}{y variable name. a continuous variable}
  \item{u}{x variable name, upper side variable}
  \item{e}{a vector of summarize function names}
  \item{ou}{order of levels of upper side x variable}
  \item{repl}{list of strings to replace after summarize. Length of list should be 2, and both should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable with one x varaible.
}
\value{
  A data.frame of summarized values. Row names are from \code{e} names. Column names are from the levels of x variable. 
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{\code{\link{tsum}}, \code{\link{tsum0}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\examples{
tsum1(CO2, "uptake", "Treatment")
tsum1(CO2, "uptake", "Treatment", 
      e=c("mean", "median", "sd", "min", "max", "length"), 
      ou=c("chilled", "nonchilled"),
      repl=list(c("median", "length"), c("med", "n")))
}
\keyword{Table Summary}
