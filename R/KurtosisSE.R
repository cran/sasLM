KurtosisSE = function(y)
{
  y = as.numeric(y)
  if (!(is.numeric(y) & is.vector(y))) return(NA)
  y = y[!is.na(y)]
  n = length(y)
  if (n < 4) return(NA) # stop("Input vector length should be larger than 3.")
  sqrt(24*n*(n - 1)^2/(n - 3)/(n - 2)/(n + 3)/(n + 5)) 
}
