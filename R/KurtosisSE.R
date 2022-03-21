KurtosisSE = function(x)
{
  x = as.numeric(x)
  if (!(is.numeric(x) & is.vector(x))) return(NA)
  x = x[!is.na(x)]
  n = length(x)
  if (n < 4) return(NA) # stop("Input vector length should be larger than 3.")
  sqrt(24*n*(n - 1)^2/(n - 3)/(n - 2)/(n + 3)/(n + 5)) 
}
