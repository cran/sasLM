\name{WhiteTest}
\alias{WhiteTest}
\title{White's Model Specification Test}
\description{
  This is shown in SAS PROC REG as the Test of First and Second Moment Specification.
}
\usage{
  WhiteTest(rx)
}
\arguments{
  \item{rx}{a result of lm}
}
\details{
  This is also called as White's general test for heteroskedasticity.
}
\value{
  Returns a direct test results by more coomplex theorem 2 , not by simpler corollary 1.
}
\examples{
  WhiteTest(lm(mpg ~ disp, mtcars))
}
\references{
  White H. A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. Econometrica 1980;48(4):817-838.
}
\author{Kyun-Seop Bae k@acr.kr}
