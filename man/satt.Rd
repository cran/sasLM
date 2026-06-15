\name{satt}
\alias{satt}
\title{Satterthwaite Approximation of Variance and Degrees of Freedom}
\description{
  Calculates the pooled variance and degrees of freedom using the Satterthwaite equation.
}
\usage{
  satt(vars, dfs, ws=c(1, 1)) 
}
\arguments{
  \item{vars}{a vector of variances}
  \item{dfs}{a vector of degrees of freedom}
  \item{ws}{a vector of weights}
}
\details{
  The input can contain more than two variances.
}
\value{
  \item{Variance}{approximated variance}
  \item{Df}{degrees of freedom}
}
\author{Kyun-Seop Bae k@acr.kr}
