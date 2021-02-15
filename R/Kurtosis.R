Kurtosis = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  d = x - mean(x)
  n*(n + 1)*sum(d^4)/((n - 1)*(n - 2)*(n - 3)*(sum(d^2)/(n - 1))^2) - 3*(n - 1)^2/((n - 2)*(n - 3))
}
