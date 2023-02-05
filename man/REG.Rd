\name{REG}
\alias{REG}
\title{Regression of Linear Least Square, similar to SAS PROC REG}
\description{
  REG is similar to SAS PROC REG.
}
\usage{
  REG(Formula, Data, conf.level=0.95, HC=FALSE, Resid=FALSE, Weights=1,
      summarize=TRUE)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{conf.level}{confidence level for the confidence limit}
  \item{HC}{heteroscedasticity related output is required such as HC0, HC3, White's first and second moment specification test}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
  \item{Weights}{weights for the each observation or residual square. This is usually the inverse of each variance.}
  \item{summarize}{If this is \code{FALSE}, \code{REG} returns just \code{lfit} result.}
}
\details{
  It performs the core function of SAS PROC REG.
}
\value{
  The result is comparable to that of SAS PROC REG.

  The first part is ANOVA table.

  The second part is measures about fitness.

  The third part is the esimates of coefficients.
  \item{Estimate}{point estimate of parameters, coefficients}
  \item{Estimable}{estimability: 1=TRUE, 0=FALSE. This appears only when at least one inestimability occurs.}
  \item{Std. Error}{standard error of the point estimate}
  \item{Lower CL}{lower confidence limit with conf.level}
  \item{Upper CL}{lower confidence limit with conf.level}
  \item{Df}{degree of freedom}
  \item{t value}{value for t distribution}
  \item{Pr(>|t|)}{probability of larger than absolute t value from t distribution with residual's degree of freedom}

  The above result is repeated using HC0 and HC3, with following White's first and second moment specification test, if HC option is specified.
  The statstic t values and their p values with HC1 and HC2 are between those of HC0 and H3.

  \item{Fitted}{Fitted value or y hat. This is returned only with Resid=TRUE option.}
  \item{Residual}{Weighted residuals. This is returned only with Resid=TRUE option.}

  If \code{summarize=FALSE}, \code{REG} returns;
  \item{coeffcients}{beta coefficients}
  \item{g2}{g2 inverse}
  \item{rank}{rank of the model matrix}
  \item{DFr}{degree of freedom for the residual}
  \item{SSE}{sum of square error}
}
\examples{
  REG(uptake ~ Plant + Type + Treatment + conc, CO2)
  REG(uptake ~ conc, CO2, HC=TRUE)
  REG(uptake ~ conc, CO2, Resid=TRUE)
  REG(uptake ~ conc, CO2, HC=TRUE, Resid=TRUE)
  REG(uptake ~ conc, CO2, summarize=FALSE)
}
\seealso{
  \code{\link{lr}}
}
\author{Kyun-Seop Bae k@acr.kr}
