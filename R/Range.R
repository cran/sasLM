Range = function(x)
{
  x = as.numeric(x)  
  if (!(is.numeric(x) & is.vector(x))) return(NA) # stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  if (length(x) == 0) return(NA)
  max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
}
