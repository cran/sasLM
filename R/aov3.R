aov3 = function(Formula, Data)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  Lx = e3(x)

  if ("Complete" %in% names(alias(Formula, Data))) {
    aNames = rownames(alias(Formula, Data)$Complete)
    if (any(colnames(Data) %in% aNames)) {
      Data[, colnames(Data) %in% aNames] = 0
      x = ModelMatrix(Formula, Data, KeepOrder=FALSE)
    }
  }

  rx = lfit(x, y)
  Tx = SS(x, rx, Lx)

  SST = as.numeric(crossprod(y - attr(x$terms, "intercept")*mean(y)))

  return(sumANOVA(rx, Tx, SST, nrow(x$X), rownames(attr(terms(x), "factors"))[1]))
}
