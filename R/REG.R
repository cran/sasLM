REG = function(Formula, Data, eps=1e-8, summarize=TRUE)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  if (summarize) {
    return(sumREG(lfit(x, y, eps=eps), x$X))
  } else {
    return(lfit(x, y, eps=eps))
  }
}
