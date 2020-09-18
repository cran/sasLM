\name{satt}
\alias{satt}
\title{Satterthwaite Approximation of Pooled Variance and Degree of Freedom}
\description{
  Calculates pooled variance and degree of freedom using Satterthwaite equation. 
}
\usage{
  satt(ws, vars, dfs) 
}
\arguments{
  \item{ws}{a vector of weights}
  \item{vars}{a vector of variances}
  \item{dfs}{a vector of degree of freedoms}
}
\details{
  The input can be more than two variances.
}
\value{
  \item{Variance}{pooled variance}
  \item{Df}{degree of freedom}
}
\author{Kyun-Seop Bae k@acr.kr}
