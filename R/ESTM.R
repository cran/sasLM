ESTM = function(L, Formula, Data, conf.level=0.95)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)

  rx = REG(Formula, Data, summarize=FALSE)
  return(est(L, x$X, rx, conf.level=conf.level))
}

