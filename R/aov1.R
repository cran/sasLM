aov1 = function(Formula, Data)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x0 = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  Lx = e1(crossprod(x0$X))
  rx = lfit(x0, y)

  x1 = x0
  fAdj = FALSE
  if ("Complete" %in% names(alias(Formula, Data))) {
    aNames = rownames(alias(Formula, Data)$Complete)
    if (any(colnames(Data) %in% aNames)) {
      Data[, colnames(Data) %in% aNames] = 0
      x1 = ModelMatrix(Formula, Data, KeepOrder=TRUE)
      Lx = e1(crossprod(x1$X))
      rx = lfit(x1, y)
      fAdj = TRUE
    }
  }

  if (fAdj) rx$coefficients = rx$coefficients*estmb(diag(NCOL(x1$X)), x1$X, rx$g2)
  Tx = SS(x1, rx, Lx)

  SST = as.numeric(crossprod(y - attr(x0$terms, "intercept")*mean(y)))

  return(sumANOVA(rx, Tx, SST, nrow(x0$X), rownames(attr(terms(x0), "factors"))[1]))
}
