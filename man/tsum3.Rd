\name{tsum3}
\alias{tsum3}
\title{Table Summary 3 independent(x) variables}
\description{
  Summarize a continuous dependent(y) variable with three independent(x) variables.   
}
\usage{
  tsum3(d, y, l, u, e=c("mean", "sd", "length"), h=NULL, ol1="", ol2="", ou="", 
        rm.dup=TRUE, repl=list(c("length"), c("n")))
}
\arguments{
  \item{d}{a data.frame or matrix with colnames}
  \item{y}{y variable name. a continuous variable}
  \item{l}{a vector of two x variable name to be shown on the left side. The length should be 2.}
  \item{u}{x variable name to be shown on the upper side}
  \item{e}{a vector of summarize function names}
  \item{h}{a list of two vectors of summarize function names for the first and second horizontal subgroups. If \code{NULL}, it becomes same to \code{e} argument.}
  \item{ol1}{order of levels of 1st left side x variable}
  \item{ol2}{order of levels of 2nd left side x variable}
  \item{ou}{order of levels of upper side x variable}
  \item{rm.dup}{if \code{TRUE}, duplicated name of levels are specified on the first occurrence only.}
  \item{repl}{list of strings to replace after summarize. Length of list should be 2, and both should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable with three x varaibles; two on the left side, the other on the upper side.
}
\value{
  A data.frame of summarized values. Column names are from the levels of \code{u}. Row names are basically from the levels of \code{l}. 
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{\code{\link{tsum}}, \code{\link{tsum0}}, \code{\link{tsum1}}, \code{\link{tsum2}}}
\examples{
tsum3(CO2, "uptake", c("Type", "Treatment"), "conc")
tsum3(CO2, "uptake", c("Type", "Treatment"), "conc", 
      e=c("mean", "median", "sd", "min", "max", "length"),
      h=list(c("mean", "sd", "length"), c("mean", "length")),
      ol2=c("chilled", "nonchilled"),
      repl=list(c("median", "length"), c("med", "n")))
}
\keyword{Table Summary}
