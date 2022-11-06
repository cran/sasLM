Median = function(y)
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  n = length(y)
  if (n == 0) return(NA)
  median(y, na.rm=TRUE)
}
