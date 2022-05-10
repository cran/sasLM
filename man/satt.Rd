\name{satt}
\alias{satt}
\title{Satterthwaite Approximation of Variance and Degree of Freedom}
\description{
  Calculates pooled variance and degree of freedom using Satterthwaite equation. 
}
\usage{
  satt(vars, dfs, ws=c(1, 1)) 
}
\arguments{
  \item{vars}{a vector of variances}
  \item{dfs}{a vector of degree of freedoms}
  \item{ws}{a vector of weights}
}
\details{
  The input can be more than two variances.
}
\value{
  \item{Variance}{approximated variance}
  \item{Df}{degree of freedom}
}
\author{Kyun-Seop Bae k@acr.kr}
