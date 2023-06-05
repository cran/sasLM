\name{pResD}
\alias{pResD}
\title{Residual Diagnostic Plot for Regression}
\description{
  Nine residual diagnostics plots.
}
\usage{
  pResD(rx, Title=NULL)
}
\arguments{
  \item{rx}{a result of lm, which can give \code{fitted}, \code{residuals}, and \code{rstandard}.}
  \item{Title}{title to be printed on the plot}
}
\details{
  SAS-style residual diagnostic plots.
}
\value{
  Nine residual diagnostic plots in a page.
}
\examples{
  pResD(lm(uptake ~ Plant + Type + Treatment + conc, CO2), "Residual Diagnostic Plot")
}
\keyword{Plot}
\author{Kyun-Seop Bae k@acr.kr}

