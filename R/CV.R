CV = function(x)
{
  x = x[!is.na(x)]
  if (length(x) == 0) return(NA)
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE) * 100
}
