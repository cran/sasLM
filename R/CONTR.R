CONTR = function(L, Formula, Data, mu = 0)
{
  rx = REG(Formula, Data, summarize=FALSE)
  if (ncol(L) != length(rx$coefficients)) stop("Check the contrast matrix!")
  return(cSS(L, rx, mu=mu))
}
