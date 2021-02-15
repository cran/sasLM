Range = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
}
