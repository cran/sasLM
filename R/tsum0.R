tsum0 = function(d, y, e=c("Mean", "SD", "N"), repl=list(c("length"), c("n")))
{
  y0 = d[!is.na(d[,y]), y]

  ne = length(e)
  Res = vector(length=ne)

  sNames = e
  for (k in 1:length(repl[[1]])) sNames[sNames == repl[[1]][k]] = repl[[2]][k]
  names(Res) = sNames

  for (k in 1:ne) Res[k] = do.call(e[k], list(y0))

  attr(Res, "y") = y
  return(Res)
}
