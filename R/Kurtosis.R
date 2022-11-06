Kurtosis = function(y)
{
  y = as.numeric(y)
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  n = length(y)
  if (n < 4) return(NA) # stop("Input vector length should be larger than 3.")
  d = y - mean(y)
  n*(n + 1)*sum(d^4)/((n - 1)*(n - 2)*(n - 3)*(sum(d^2)/(n - 1))^2) - 3*(n - 1)^2/((n - 2)*(n - 3))
}
