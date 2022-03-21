Kurtosis = function(x)
{
  x = as.numeric(x)
  if (!(is.numeric(x) & is.vector(x))) return(NA) # stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 4) return(NA) # stop("Input vector length should be larger than 3.")
  d = x - mean(x)
  n*(n + 1)*sum(d^4)/((n - 1)*(n - 2)*(n - 3)*(sum(d^2)/(n - 1))^2) - 3*(n - 1)^2/((n - 2)*(n - 3))
}
