QuartileRange = function(x, Type=2) 
{
  x = as.numeric(x)  
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")   
  IQR(x, na.rm=TRUE, type=Type)
}
