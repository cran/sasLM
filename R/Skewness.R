Skewness = function(y)
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  n = length(y)
  if (n < 3) return(NA) # stop("Input vector length should be larger than 2.")
  d = y - mean(y, na.rm=TRUE)
  n*sqrt(n - 1)*(sum(d^3)/((n - 2)*sum(d^2)^(3/2)))
}
