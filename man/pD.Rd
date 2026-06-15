\name{pD}
\alias{pD}
\title{Diagnostic Plot for Regression}
\description{
  Four standard diagnostic plots for regression.
}
\usage{
  pD(rx, Title=NULL)
}
\arguments{
  \item{rx}{a result of \code{lm}, which can give \code{fitted}, \code{residuals}, and \code{rstandard}.}
  \item{Title}{title to be printed on the plot}
}
\details{
  The most frequently used diagnostic plots are 'observed vs. fitted', 'standardized residual vs. fitted', 'distribution plot of standardized residuals', and 'Q-Q plot of standardized residuals'.
}
\value{
  Four diagnostic plots on a page.
}
\examples{
  pD(lm(uptake ~ Plant + Type + Treatment + conc, CO2), "Diagnostic Plot")
}
\keyword{Plot}
\author{Kyun-Seop Bae k@acr.kr}
