geoCV = function(x)
{
  x = x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA)
  lx = log(x)
  vlx = var(lx, na.rm=T)
  return(sqrt(exp(vlx) - 1)*100)
}
