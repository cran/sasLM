geoMean = function(x)
{
  x = x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA)
  lx = log(x)
  return(exp(mean(lx, na.rm=T)))
}
