SkewnessSE = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 3) stop("Input vector length should be larger than 2.")
  sqrt(6*n*(n - 1)/(n - 2)/(n + 1)/(n + 3)) 
}
