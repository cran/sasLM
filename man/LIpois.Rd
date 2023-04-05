\name{LIpois}
\alias{LIpois}
\title{Likelihood Interval of the Mean assuming Poisson Distribution}
\description{
  Likelihood interval of lambda assuming Poisson distribution.
}
\usage{
  LIpois(x, k, conf.level=0.95, eps=1e-8)
}
\arguments{
  \item{x}{mean or lambda, the count in a time unit.}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
  \item{eps}{Values less than eps are considered as 0.}
}
\details{
  It calculates likelihood interval of mean(lambda) assuming Poisson distribution. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{point estimation for the lambda}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  LIpois(4, k=1/0.15) # Pawitan
  LIpois(4, k=exp(2)) # Edwards
  LIpois(4, k=8)      # Rhode
  LIpois(4)           # Bae
  LIpois(4, k=15)     # Fisher
  poisson.test(4)
  LIpois(4, k=32)     # 0.7454614 11.7893612
}
\author{Kyun-Seop Bae k@acr.kr}
