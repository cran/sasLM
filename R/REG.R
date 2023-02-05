REG = function(Formula, Data, conf.level=0.95, HC=FALSE, Resid=FALSE, Weights=1, summarize=TRUE)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  mf0 = model.frame(Formula, Data)
  y = mf0[, 1]
  n = length(y)
  yName = colnames(mf0)[1]
  if (!is.numeric(y)) stop("Response variable should be numeric!")

  x0 = ModelMatrix(Formula, Data, KeepOrder=TRUE)

  if ("Complete" %in% names(alias(Formula, Data))) {
    warning("Complete aliased variable(s) exist(s)!")
    Data[, colnames(Data) %in% rownames(alias(Formula, Data)$Complete)] = 0  # THIS IS NECESSARY
  }

  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)

  if (length(Weights) == 1 & Weights[1] == 1) Weights = rep(1, n)
  rw = sqrt(Weights) ##
  x$X = rw*x$X       ##
  y = rw*y          ##
  r0 = lfit(x, y)

  if (summarize) {
#    p = length(r0$rank)

## ANOVA
    Lx = e1(crossprod(x$X))
    Tx = SS(x, r0, Lx)
    SST = sum(Tx[Tx[, "Df"] > 0 ,"Sum Sq"]) + r0$SSE ##
    ANOVA = sumANOVA(r0, T1=NULL, SST, n, yName)

## FITNESS
    fIntercept = attr(x$terms, "intercept")
    MeanY = sum(Weights*mf0[, 1], na.rm=T)/sum(Weights) # original scale

    Rsq = ANOVA["MODEL", "Sum Sq"]/ANOVA[3, "Sum Sq"] # Corrected total or uncorrected total
    RMSE = sqrt(ANOVA["RESIDUALS", "Mean Sq"])
    CV = 100*RMSE/MeanY

    yhat = as.numeric(x$X %*% r0$coefficients)
    Residual = y - yhat
    we2 = (Residual)^2 # we2 = weighted residual^2
    iXpX = r0$g2
    hii = diag(x$X %*% iXpX %*% t(x$X))
    PRESS = sum(we2/(1 - hii)^2)
    R2pred = 1 - PRESS/SST

    if (r0$DFr > 0) {
      Rsq.adj = 1 - (1 - Rsq)*(n - fIntercept)/r0$DFr
      Fit = data.frame(RMSE, MeanY, CV, Rsq, Rsq.adj, PRESS, R2pred)
      colnames(Fit) = c("Root MSE", paste(yName, "Mean"), "Coef Var", "R-square", "Adj R-sq", "PRESS", "R2pred")
    } else {
      Fit = data.frame(RMSE, MeanY, CV, Rsq, PRESS, R2pred)
      colnames(Fit) = c("Root MSE", paste(yName, "Mean"), "Coef Var", "R-square", "PRESS", "R2pred")
    }
    rownames(Fit) = ""

    Res0 = sumREG(r0, rw*x0$X)
#    pIDX = (!is.na(Res0[, "Std. Error"]) & (Res0[, "Std. Error"] > 0))
#    ConfBound = qt(0.5 + conf.level/2, Res0[pIDX, "Df"])*Res0[pIDX, "Std. Error"]
#    LCL0 = rep(NA_real_, NROW(Res0))
#    UCL0 = LCL0
#    LCL0[pIDX] = Res0[pIDX, "Estimate"] - ConfBound
#    UCL0[pIDX] = Res0[pIDX, "Estimate"] + ConfBound
    ConfBound = qt(0.5 + conf.level/2, Res0[, "Df"])*Res0[, "Std. Error"]
    LCL0 = Res0[, "Estimate"] - ConfBound
    UCL0 = Res0[, "Estimate"] + ConfBound
    Res0 = cbind(Res0, 'Lower CL'=LCL0, 'Upper CL'=UCL0)
    Res0 = Res0[, c("Estimate", "Std. Error", "Df", "Lower CL", "Upper CL", "t value", "Pr(>|t|)"), drop=F]
    class(Res0) = "anova"

    Result = list(ANOVA=ANOVA, Fitness=Fit, Coefficients=Res0)
    iNext = 4

## HC
    if (HC) {
      Res1 = Res0
#    Res2 = Res0
#    Res3 = Res0
      Res4 = Res0

      HC0 = iXpX %*% (t(x$X) %*% diag(we2, nrow=n) %*% x$X) %*% iXpX
#    HC1 = n/(n - p)*HC0
#    HC2 = iXpX %*% (t(x$X) %*% diag(we2/(1 - hii)) %*% x$X) %*% iXpX
      HC3 = iXpX %*% (t(x$X) %*% diag(we2/(1 - hii)^2, nrow=n) %*% x$X) %*% iXpX

      Res1[, "Std. Error"] = sqrt(diag(HC0))
#    Res2[, "Std. Error"] = sqrt(diag(HC1))
#    Res3[, "Std. Error"] = sqrt(diag(HC2))
      Res4[, "Std. Error"] = sqrt(diag(HC3))

      Res1[, "t value"] = Res1[, "Estimate"]/Res1[, "Std. Error"]
#    Res2[, "t value"] = Res2[, "Estimate"]/Res2[, "Std. Error"]
#    Res3[, "t value"] = Res3[, "Estimate"]/Res3[, "Std. Error"]
      Res4[, "t value"] = Res4[, "Estimate"]/Res4[, "Std. Error"]

      Res1[, "Pr(>|t|)"] = 2*(1 - pt(abs(Res1[, "t value"]), Res1[, "Df"]))
#    Res2[, "Pr(>|t|)"] = 2*(1 - pt(abs(Res2[, "t value"]), Res2[, "Df"]))
#    Res3[, "Pr(>|t|)"] = 2*(1 - pt(abs(Res3[, "t value"]), Res3[, "Df"]))
      Res4[, "Pr(>|t|)"] = 2*(1 - pt(abs(Res4[, "t value"]), Res4[, "Df"]))

      pIDX1 = (!is.na(Res1[, "Std. Error"]) & (Res1[, "Std. Error"] > 0))
      ConfBound1 = qt(0.5 + conf.level/2, Res1[pIDX1, "Df"])*Res1[pIDX1, "Std. Error"]
      Res1[pIDX1, "Lower CL"] = Res1[pIDX1, "Estimate"] - ConfBound1
      Res1[pIDX1, "Upper CL"] = Res1[pIDX1, "Estimate"] + ConfBound1

      pIDX4 = (!is.na(Res4[, "Std. Error"]) & (Res4[, "Std. Error"] > 0))
      ConfBound4 = qt(0.5 + conf.level/2, Res4[pIDX4, "Df"])*Res4[pIDX4, "Std. Error"]
      Res4[pIDX4, "Lower CL"] = Res4[pIDX4, "Estimate"] - ConfBound4
      Res4[pIDX4, "Upper CL"] = Res4[pIDX4, "Estimate"] + ConfBound4

      Res1[is.nan(Res1)] = NA
#    Res2[is.nan(Res2)] = NA
#    Res3[is.nan(Res3)] = NA
      Res4[is.nan(Res4)] = NA

## WHITE
      Res5 = WhiteH(x$X, we2)

#    Result = list(Default=Res0, HC0=Res1, HC1=Res2, HC2=Res3, HC3=Res4, WhiteTest=Res5)
      Result[[iNext]] = Res1
      Result[[iNext + 1]] = Res4
      Result[[iNext + 2]] = Res5
      names(Result)[iNext:(iNext + 2)] = c("HC0", "HC3", "White Test")
      iNext = iNext + 3
    }

    if (Resid) {
      Result[[iNext]] = yhat
      names(Result)[iNext] = "Fitted"

      Result[[iNext + 1]] = Residual
      names(Result)[iNext + 1] = "Residual"
    }

    return(Result)
  } else {
    return(r0)
  }
}
