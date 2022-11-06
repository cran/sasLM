N = function(y)
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  length(y)
}
