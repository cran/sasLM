\name{tsum2}
\alias{tsum2}
\title{Table Summary 2 independent(x) variables}
\description{
  Summarize a continuous dependent(y) variable with two independent(x) variables.   
}
\usage{
  tsum2(d, y, l, u, e=c("Mean", "SD", "N"), h=NULL, ol="", ou="", rm.dup=TRUE, 
        repl=list(c("length"), c("n")))
}
\arguments{
  \item{d}{a data.frame or matrix with colnames}
  \item{y}{y variable name. a continuous variable}
  \item{l}{x variable name to be shown on the left side}
  \item{u}{x variable name to be shown on the upper side}
  \item{e}{a vector of summarize function names}
  \item{h}{a vector of summarize function names for the horizontal subgroup. If \code{NULL}, it becomes same to \code{e} argument.}
  \item{ol}{order of levels of left side x variable}
  \item{ou}{order of levels of upper side x variable}
  \item{rm.dup}{if \code{TRUE}, duplicated name of levels are specified on the first occurrence only.}
  \item{repl}{list of strings to replace after summarize. Length of list should be 2, and both should have the same length.}
}
\details{
  A convenient summarization function for a continuous variable with two x varaibles; one on the left side, the other on the upper side.
}
\value{
  A data.frame of summarized values. Column names are from the levels of \code{u}. Row names are basically from the levels of \code{l}. 
}
\author{Kyun-Seop Bae k@acr.kr}
\keyword{Table Summary}
\seealso{\code{\link{tsum}}, \code{\link{tsum0}}, \code{\link{tsum1}}, \code{\link{tsum3}}}
\examples{
  tsum2(CO2, "uptake", "Type", "Treatment")
  tsum2(CO2, "uptake", "Type", "conc")
  tsum2(CO2, "uptake", "Type", "Treatment", 
        e=c("mean", "median", "sd", "min", "max", "length"), 
        ou=c("chilled", "nonchilled"),
        repl=list(c("median", "length"), c("med", "n")))
}
