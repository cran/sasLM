geoMean = function(y)
{
  y = y[!is.na(y) & y > 0]
  if (length(y) == 0) return(NA)
  ly = log(y)
  return(exp(mean(ly, na.rm=T)))
}
