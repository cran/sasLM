satt = function(ws, vars, dfs) 
{
  Variance = sum(ws*vars)
  Df = Variance^2/sum((ws*vars)^2/dfs)
  Result = list(Variance=Variance, Df=Df)
  return(Result)
}