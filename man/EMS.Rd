\name{EMS}
\alias{EMS}
\title{Expected Mean Square Formula}
\description{
  Calculates a formula table for expected mean square of the given contrast. The default is for Type III SS.
}
\usage{
  EMS(Formula, Data, Type=3, eps=1e-8)
}
\arguments{
  \item{Formula}{a conventional formula for a linear model}
  \item{Data}{a \code{data.frame} to be analyzed}
  \item{Type}{type of sum of squares. The default is 3. Type 4 is not supported yet.}
  \item{eps}{Less than this value is considered as zero.}  
}
\details{
  This is necessary for further hypothesis tests of nesting factors.
}
\value{
  A coefficient matrix for Type III expected mean square
}
\author{Kyun-Seop Bae k@acr.kr}
\examples{
  f1 = log(CMAX) ~ SEQ/SUBJ + PRD + TRT
  EMS(f1, BEdata)
  EMS(f1, BEdata, Type=1)
  EMS(f1, BEdata, Type=2)
}
