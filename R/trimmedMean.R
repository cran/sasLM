trimmedMean = function(y, Trim=0.05) 
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) stop("Input should be a numeric vector!")  
  mean(y, trim=Trim, na.rm=TRUE)
}

