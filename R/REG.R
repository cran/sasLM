REG = function(Formula, Data, NOINT=FALSE, eps=1e-8)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE, NOINT=NOINT)
  return(sumREG(lfit(x, y, eps=eps), x$X))
}
