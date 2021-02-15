Mean = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 1) stop("The length of the input vector should be at least 1.")
  mean(x)
}
