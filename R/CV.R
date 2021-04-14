CV = function(x)
{
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE) * 100
}
