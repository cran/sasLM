geoCV = function(y)
{
  y = y[!is.na(y) & y > 0]
  if (length(y) == 0) return(NA)
  ly = log(y)
  vly = var(ly, na.rm=T)
  return(sqrt(exp(vly) - 1)*100)
}
