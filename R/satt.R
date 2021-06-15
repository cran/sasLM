satt = function(vars, dfs, ws=c(1, 1)) 
{
  if (length(vars) != length(dfs)) stop("Lengths of vars and dfs are different!")
  if (length(ws) != length(vars)) ws = rep(1, length(vars))
  ws = ws/sum(ws)
  Variance = sum(ws*vars)
  Df = Variance^2/sum((ws*vars)^2/dfs)
  Result = list(Variance=Variance, Df=Df)
  return(Result)
}
