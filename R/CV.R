CV = function(y)
{
  y = y[!is.na(y)]
  if (length(y) == 0) return(NA)
  sd(y, na.rm=TRUE)/mean(y, na.rm=TRUE)*100
}
