SkewnessSE = function(y)
{
  y = as.numeric(y)
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  n = length(y)
  if (n < 3) return(NA) # stop("Input vector length should be larger than 2.")
  sqrt(6*n*(n - 1)/(n - 2)/(n + 1)/(n + 3)) 
}
