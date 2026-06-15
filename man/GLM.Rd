\name{GLM}
\alias{GLM}
\title{General Linear Model similar to SAS PROC GLM}
\description{
  GLM is the main function of this package.
}
\usage{
  GLM(Formula, Data, BETA=FALSE, EMEAN=FALSE, Resid=FALSE, conf.level=0.95,
      Weights=1)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{BETA}{if \code{TRUE}, coefficients (parameters) of \code{REG} will be returned. This is equivalent to the SOLUTION option of SAS PROC GLM}
  \item{EMEAN}{if \code{TRUE}, least square means (or expected means) will be returned. This is equivalent to the LSMEANS statement of SAS PROC GLM}
  \item{Resid}{if \code{TRUE}, fitted values (y hat) and residuals will be returned}
  \item{conf.level}{confidence level for the confidence limit of the least square mean}
  \item{Weights}{weights for the weighted least squares. This should be a scalar or a vector of the same length as the number of rows of \code{Data}. Observations with nonpositive weights are excluded from the analysis as SAS does.}
}
\details{
  It performs the core function of SAS PROC GLM. Least square means for the interaction term of three variables are not supported yet.
}
\value{
  The result is comparable to that of SAS PROC GLM.
  \item{ANOVA}{ANOVA table for the model}
  \item{Fitness}{Some measures of goodness of fit such as R-square and CV}
  \item{Type I}{Type I sum of squares table}
  \item{Type II}{Type II sum of squares table}
  \item{Type III}{Type III sum of squares table}
  \item{Parameter}{Parameter table with standard error, t value, p value. \code{TRUE} is 1, and \code{FALSE} is 0 in the \code{Estimable} column. This is returned only with the BETA=TRUE option.}
  \item{Expected Mean}{Least square (or expected) mean table with confidence limits. This is returned only with the EMEAN=TRUE option.}
  \item{Fitted}{Fitted values or y hat in the original scale, as SAS OUTPUT P= does, even with Weights. This is returned only with the Resid=TRUE option.}
  \item{Residual}{Residuals in the original scale, as SAS OUTPUT R= does, even with Weights. This is returned only with the Resid=TRUE option.}
}
\examples{
GLM(uptake ~ Type*Treatment + conc, CO2[-1,]) # Making data unbalanced
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], EMEAN=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], Resid=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE, EMEAN=TRUE)
GLM(uptake ~ Type*Treatment + conc, CO2[-1,], BETA=TRUE, EMEAN=TRUE, Resid=TRUE)
}
\seealso{
  \code{\link{REG}}, \code{\link{aov1}}, \code{\link{aov2}}, \code{\link{aov3}}, \code{\link{LSM}}, \code{\link{PDIFF}}
}
\author{Kyun-Seop Bae k@acr.kr}
