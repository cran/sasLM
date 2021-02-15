\name{bk}
\alias{bk}
\title{Beautify the output of knitr::kable}
\description{
  Trailing zeros after integer is somwwhat annoying. This removes those in the vector of strings.   
}
\usage{
  bk(ktab, rpltag=c("n", "N"), dig=10)
}
\arguments{
  \item{ktab}{an output of \code{knitr::kable}}
  \item{rpltag}{tag string of replacement rows. This is usually "n" which means the sample count.}
  \item{dig}{maximum digits of decimals in the \code{kable} output}
}
\details{
  This is convenient if used with tsum0, tsum1, tsum2, tsum3, This requires \code{knitr::kable}.
}
\value{
  A new processed vector of strings. The class is still \code{knitr_kable}.
}
\author{Kyun-Seop Bae k@acr.kr}
\seealso{\code{\link{tsum0}}, \code{\link{tsum1}}, \code{\link{tsum2}}, \code{\link{tsum3}}}
\examples{
## OUTPUT example
# t0 = tsum0(CO2, "uptake", c("mean", "median", "sd", "length", "min", "max"))
# bk(kable(t0)) # requires knitr package
#
# |       |        x|
# |:------|--------:|
# |mean   | 27.21310|
# |median | 28.30000|
# |sd     | 10.81441|
# |n      | 84      |
# |min    |  7.70000|
# |max    | 45.50000|

# t1 = tsum(uptake ~ Treatment, CO2, 
#           e=c("mean", "median", "sd", "min", "max", "length"), 
#           ou=c("chilled", "nonchilled"),
#           repl=list(c("median", "length"), c("med", "N")))
# 
# bk(kable(t1, digits=3)) # requires knitr package
# 
# |     | chilled| nonchilled| Combined|
# |:----|-------:|----------:|--------:|
# |mean |  23.783|     30.643|   27.213|
# |med  |  19.700|     31.300|   28.300|
# |sd   |  10.884|      9.705|   10.814|
# |min  |   7.700|     10.600|    7.700|
# |max  |  42.400|     45.500|   45.500|
# |N    |  42    |     42    |   84    |
}


