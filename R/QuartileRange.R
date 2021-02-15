QuartileRange = function(x, Type=6) 
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")   
  IQR(x, na.rm=TRUE, type=Type)
}
