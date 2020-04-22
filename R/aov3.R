aov3 = function(Formula, Data)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  r1 = lfit(x, y)
  T1 = SS(x, r1, e3(Formula, Data))
  return(sumANOVA(r1, T1, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x),"factors"))[1]))
}
