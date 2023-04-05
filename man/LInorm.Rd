\name{LInorm}
\alias{LInorm}
\title{Likelihood Interval of mean and sd assuming Norml Distribution}
\description{
  Likelihood interval of mean and sd assuming normal distribution. This is estimated likelihood interval, not profile likelihood interval.
}
\usage{
  LInorm(x, k, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{x}{a vector of observation}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
  \item{eps}{Values less than eps are considered as 0.}
}
\details{
  It calculates likelihood interval of mean and sd assuming Poisson distribution in one group. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{point estimation for the proportion}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  x = c(-5.3, -4.5, -1.0, -0.7, 3.7, 3.9, 4.2, 5.5, 6.8, 7.4, 9.3)
  LInorm(x, k=1/0.15) # Pawitan Ex10-9 p289
  LInorm(x)
}
\author{Kyun-Seop Bae k@acr.kr}
