Skewness = function(x)
{
  x = as.numeric(x)  
  if (!(is.numeric(x) & is.vector(x))) return(NA) # stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 3) return(NA) # stop("Input vector length should be larger than 2.")
  d = x - mean(x, na.rm=TRUE)
  n*sqrt(n - 1)*(sum(d^3)/((n - 2)*sum(d^2)^(3/2)))
}
