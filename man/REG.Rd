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
  \item{HC}{heteroscedasticity-related output is required, such as HC0, HC3, and White's first and second moment specification test}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
  \item{Weights}{weights for each observation, usually the inverse of each variance. This should be a scalar or a vector of the same length as the number of rows of \code{Data}. Observations with nonpositive weights are excluded from the analysis, as SAS does.}
  \item{summarize}{If this is \code{FALSE}, \code{REG} returns just the \code{lfit} result.}
}
\details{
  It performs the core function of SAS PROC REG.
}
\value{
  The result is comparable to that of SAS PROC REG.

  The first part is the ANOVA table.

  The second part is measures of fitness.

  The third part is the estimates of coefficients.
  \item{Estimate}{point estimate of parameters, coefficients}
  \item{Estimable}{estimability: 1=TRUE, 0=FALSE. This appears only when at least one inestimability occurs.}
  \item{Std. Error}{standard error of the point estimate}
  \item{Lower CL}{lower confidence limit with conf.level}
  \item{Upper CL}{upper confidence limit with conf.level}
  \item{Df}{degrees of freedom}
  \item{t value}{value for the t distribution}
  \item{Pr(>|t|)}{probability of a larger absolute t value from the t distribution with the residual's degrees of freedom}

  The above result is repeated using HC0 and HC3, followed by White's first and second moment specification test, if the HC option is specified.
  The t values and their p values with HC1 and HC2 are between those of HC0 and HC3.

  \item{Fitted}{Fitted value or y hat in the original scale as SAS OUTPUT P= does, even with Weights. This is returned only with the Resid=TRUE option.}
  \item{Residual}{Residuals in the original scale as SAS OUTPUT R= does, even with Weights. This is returned only with the Resid=TRUE option.}

  If \code{summarize=FALSE}, \code{REG} returns;
  \item{coefficients}{beta coefficients}
  \item{g2}{g2 inverse}
  \item{rank}{rank of the model matrix}
  \item{DFr}{degrees of freedom for the residual}
  \item{SSE}{sum of squared errors}
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
