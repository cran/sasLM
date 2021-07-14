SD = function(x)
{
  if (!(is.numeric(x) & is.vector(x))) return(NA) # stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 2) return(NA)
  sd(x, na.rm=TRUE)
}
