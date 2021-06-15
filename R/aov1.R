aov1 = function(Formula, Data, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  r1 = lfit(x, y, eps=eps)
  T1 = SS(x, r1, e1(Formula, Data, eps=eps))
  if ("(Intercept)" %in% colnames(x$X)) {
    SST = crossprod(y - mean(y))
  } else {
    SST = crossprod(y)
  }
  return(sumANOVA(r1, T1, SST, nrow(x$X), rownames(attr(terms(x),"factors"))[1]))
}
