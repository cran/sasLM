Range = function(y)
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  if (length(y) == 0) return(NA)
  max(y, na.rm=TRUE) - min(y, na.rm=TRUE)
}
