QuartileRange = function(y, Type=2) 
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) stop("Input should be a numeric vector!")   
  IQR(y, na.rm=TRUE, type=Type)
}
