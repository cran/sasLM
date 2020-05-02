\name{ANOVA}
\alias{ANOVA}
\title{Analysis of Variance similar to SAS PROC ANOVA}
\description{
  Analysis of variance with type I, II, and III sum of squares.
}
\usage{
ANOVA(Formula, Data, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model.}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{eps}{Less than this value is considered as zero.}
}
\details{
  It performs the core function of SAS PROC ANOVA.
}
\value{
  The result is comparable to that of SAS PROC ANOVA.
  \item{ANOVA}{ANOVA table for the model}
  \item{Type I}{Type I sum of square table}
  \item{Type II}{Type II sum of square table}
  \item{Type III}{Type III sum of square table}
}
\examples{
  ANOVA(uptake ~ Plant + Type + Treatment + conc, CO2)
}
\author{Kyun-Seop Bae k@acr.kr}
