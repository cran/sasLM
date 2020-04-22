af = function(DataFrame, Cols)
{
  for (i in Cols) DataFrame[,i] = as.factor(DataFrame[,i])
  return(DataFrame)
}
