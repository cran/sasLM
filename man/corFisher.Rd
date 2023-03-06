\name{corFisher}
\alias{corFisher}
\title{Correlation test by Fisher's Z transformation}
\description{
  Testing correlation between two numeric vectors by Fisher's Z transformation
}
\usage{
  corFisher(x, y, conf.level=0.95, rho=0) 
}
\arguments{
  \item{x}{the first input numeric vector}
  \item{y}{the second input numeric vector}
  \item{conf.level}{confidence level}
  \item{rho}{population correlation rho under null hypothesis}  
}
\details{
  This accepts only two numeric vectors.
}
\value{
  \item{N}{sample size, length of input vectors}
  \item{r}{sample correlation}
  \item{Fisher.z}{Fisher's z}
  \item{bias}{bias to correct}
  \item{rho.hat}{point estimate of population rho}
  \item{conf.level}{confidence level for the confidence interval}
  \item{lower}{lower limit of confidence interval}
  \item{upper}{upper limit of confidence interval}
  \item{rho0}{population correlation rho under null hypothesis}
  \item{p.value}{p value under the null hypothesis}
}
\author{Kyun-Seop Bae k@acr.kr}
\references{Fisher RA. Statistical Methods for Research Workers. 14e. 1973}
\examples{
  corFisher(mtcars$disp, mtcars$hp, rho=0.6)
}
