\name{TTEST}
\alias{TTEST}
\title{Independent two groups t-test comparable to PROC TTEST}
\description{
  This is comparable to SAS PROC TTEST.
}
\usage{
  TTEST(x, y, conf.level=0.95)
}
\arguments{
  \item{x}{a vector of data from the first (test, active, experimental) group}
  \item{y}{a vector of data from the second (reference, control, placebo) group}
  \item{conf.level}{confidence level}
}
\details{
  Caution on choosing the row to use in the output.
}
\value{
  The output format is comparable to SAS PROC TTEST.
}
\examples{
  TTEST(mtcars[mtcars$am==1, "mpg"], mtcars[mtcars$am==0, "mpg"])
}
\seealso{\code{\link{mtest}}, \code{\link{tmtest}}, \code{\link{ztest}}}
\keyword{t-test}
\author{Kyun-Seop Bae k@acr.kr}
