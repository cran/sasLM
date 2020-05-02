aov1 = function(Formula, Data, eps=1e-8)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  r1 = lfit(x, y, eps=eps)
  T1 = SS(x, r1, e1(Formula, Data, eps=eps))
  return(sumANOVA(r1, T1, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x),"factors"))[1]))
}
