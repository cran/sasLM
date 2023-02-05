aov2 = function(Formula, Data, BETA=FALSE, Resid=FALSE)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x0 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  Lx = e2(x0)
  rx = lfit(x0, y)
  vESTM = estmb(diag(NCOL(x0$X)), x0$X, rx$g2)

  x1 = x0
  fAdj = FALSE
  if ("Complete" %in% names(alias(Formula, Data))) {
    aNames = rownames(alias(Formula, Data)$Complete)
    if (any(colnames(Data) %in% aNames)) {
      Data[, colnames(Data) %in% aNames] = 0
      x1 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
      rx = lfit(x1, y)
      fAdj = TRUE
    }
  }

  if (fAdj) rx$coefficients = rx$coefficients*vESTM
  Tx = SS(x1, rx, Lx)
  SST = as.numeric(crossprod(y - attr(x0$terms, "intercept")*mean(y)))
  Res0 = sumANOVA(rx, Tx, SST, nrow(x0$X), rownames(attr(terms(x0), "factors"))[1])

  if (!BETA & !Resid) {
    Result = Res0
  } else {
    Result = list(ANOVA=Res0)
    iNext = 2

    if (BETA) {
      Result[[iNext]] = sumREG(rx, x1$X)
      names(Result)[iNext] = "Parameter"
      iNext = iNext + 1
    }

    if (Resid) {
      yhat = as.numeric(x1$X %*% rx$coefficients)
      Result[[iNext]] = yhat
      names(Result)[iNext] = "Fitted"

      Result[[iNext + 1]] = y - yhat
      names(Result)[iNext + 1] = "Residual"
    }
  }

  return(Result)
}
