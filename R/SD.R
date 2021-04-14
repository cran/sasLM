SD = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 2) stop("The length of the input vector should be larger than 1.")
  sd(x)
}