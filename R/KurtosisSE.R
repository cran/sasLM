KurtosisSE = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 4) stop("Input vector length should be larger than 3.")
  sqrt(24*n*(n - 1)^2/(n - 3)/(n - 2)/(n + 3)/(n + 5)) 
}
