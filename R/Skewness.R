Skewness = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  d = x - mean(x)
  n*sqrt(n - 1)*(sum(d^3)/((n - 2)*sum(d^2)^(3/2)))
}
