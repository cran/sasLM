trimmedMean = function(x, Trim=0.05) 
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")  
  mean(x, trim=Trim, na.rm=TRUE)
}

