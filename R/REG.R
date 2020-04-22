REG = function(Formula, Data, NOINT=FALSE)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE, NOINT=NOINT)
  return(sumREG(lfit(x, y), x$X))
}
