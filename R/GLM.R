GLM = function(Formula, Data, BETA=FALSE, EMEAN=FALSE, Resid=FALSE, conf.level=0.95, Weights=1)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  mf0 = model.frame(Formula, Data)
  nObs = NROW(mf0)
  srwi = sqrt(Weights) ##
  y = srwi*mf0[, 1] ##
  yName = names(mf0)[1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x0a = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  x0b = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  x0a$X = srwi*x0a$X ##
  x0b$X = srwi*x0b$X ##
  r1 = lfit(x0a, y)
  r2 = lfit(x0b, y)
  r3 = r2
  rex = ex(x0a, x0b, r2$g2)
  nPara = NCOL(x0a$X)
  vESTM = estmb(diag(nPara), x0a$X, r1$g2)

  x1 = x0a
  L1 = rex$e1
  fAdj = FALSE
  if ("Complete" %in% names(alias(Formula, Data))) {
    aNames = rownames(alias(Formula, Data)$Complete)
    if (any(colnames(Data) %in% aNames)) {
      Data[, colnames(Data) %in% aNames] = 0
      x1 = ModelMatrix(Formula, Data, KeepOrder=TRUE)
      x2 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
      x1$X = srwi*x1$X ##
      x2$X = srwi*x2$X ##
      r1 = lfit(x1, y)
      r2 = lfit(x2, y)
      L1 = e1(crossprod(x1$X))
      fAdj = TRUE
    }
  }

  r1x = r1
  if (fAdj) r1x$coefficients = r1x$coefficients*estmb(diag(nPara), x1$X, r1$g2)
  T1 = SS(x0a, r1x, L1)
  class(T1) = "anova"

  if (fAdj) r2$coefficients = r2$coefficients*vESTM
  T2 = SS(x0b, r2, rex$e2)[rownames(T1), , drop=FALSE]
  class(T2) = "anova"

  T3 = SS(x0b, r3, rex$e3)[rownames(T1), , drop=FALSE]
  if (sum(T3[, "Df"], na.rm=TRUE) != sum(T1[, "Df"], na.rm=TRUE)) {
    attr(T3, "heading") = "CAUTION: Singularity Exists !"
  }
  class(T3) = "anova"

  fIntercept = attr(x0a$terms, "intercept")
  if (length(Weights) == 1 & Weights[1] == 1) { MeanY = mean(mf0[, 1], na.rm=T)
  } else { MeanY = sum(Weights*mf0[, 1], na.rm=T)/sum(Weights) } ## weighted mean
  SST = sum(T1[T1[, "Df"] > 0, "Sum Sq"]) + r1x$SSE ##

  ANOVA = sumANOVA(r1, T1=NULL, SST, nObs, rownames(attr(terms(x0a), "factors"))[1])

  Rsq = ANOVA["MODEL", "Sum Sq"]/ANOVA[3, "Sum Sq"] # Corrected total or uncorrected total
  RMSE = sqrt(ANOVA["RESIDUALS", "Mean Sq"])
  CV = 100*RMSE/MeanY

  if (r1$DFr > 0) {
    Rsq.adj = 1 - (1 - Rsq)*(nObs - fIntercept)/r1$DFr
    Fit = data.frame(RMSE, MeanY, CV, Rsq, Rsq.adj)
    colnames(Fit) = c("Root MSE", paste(yName, "Mean"), "Coef Var", "R-square", "Adj R-sq")
  } else {
    Fit = data.frame(RMSE, MeanY, CV, Rsq)
    colnames(Fit) = c("Root MSE", paste(yName, "Mean"), "Coef Var", "R-square")
  }
  rownames(Fit) = ""

  Result = list(ANOVA=ANOVA, Fitness=Fit, 'Type I'=T1, 'Type II'=T2, 'Type III'=T3)
  iNext = 6

  if (BETA) {
    Result[[iNext]] = sumREG(r1, x0a$X)
    names(Result)[iNext] = "Parameter"
    iNext = iNext + 1
  }

  if (EMEAN) {
    Result[[iNext]] = lsm0(x2, r2, Formula, Data, conf.level=conf.level)
    names(Result)[iNext] = "Expected Mean"
    iNext = iNext + 1
  }

  if (Resid) {
    yhat = as.numeric(x1$X %*% r1$coefficients)
    Result[[iNext]] = yhat
    names(Result)[iNext] = "Fitted"

    Result[[iNext + 1]] = y - yhat
    names(Result)[iNext + 1] = "Residual"
  }

  return(Result)
}
