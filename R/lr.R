lr = function(Formula, Data, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  if ("Complete" %in% names(alias(Formula, Data))) {
    warning("Complete aliased variable(s) exist(s)!")
    Data[, rownames(alias(Formula, Data)$Complete)] = 0
  }

  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x = ModelMatrix(Formula, Data)
  nc = ncol(x$X)
  XpX = crossprod(x$X)
  XpY = crossprod(x$X, y)
  aXpX = rbind(cbind(XpX, XpY), cbind(t(XpY), crossprod(y)))
  ag2 = G2SWEEP(aXpX, Augmented=TRUE, eps=eps)
  b = ag2[1:nc, (nc + 1)]
  iXpX = ag2[1:nc, 1:nc, drop=FALSE]
  nr = nrow(x$X)
  np = attr(ag2, "rank")
  DFr = nr - np
  SSE = max(0, ag2[(nc + 1), (nc + 1)])
  fIntercept = attr(x$terms, "intercept")
  SST = as.numeric(crossprod(y - fIntercept*mean(y)))

  if (DFr > 0) {
    MSE = SSE/DFr
    bVar = iXpX %*% XpX %*% t(iXpX) * MSE
    bVar[abs(bVar) < eps] = NA_real_
    bSE = sqrt(diag(bVar))
    Tval = b/bSE
    Pval = 2*(1 - pt(abs(Tval), DFr))
  } else {
    MSE = NA
    bSE = NA
    Tval = NA
    Pval = NA
  }

  Parameter = cbind(b, bSE, Tval, Pval)
  colnames(Parameter) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(Parameter) = colnames(x$X)

  Res = list()
  Res$call = match.call()
  Res$terms = x$terms
  Res$residuals = as.vector(y - x$X %*% b)

  coef1 = Parameter
  coef1[is.na(bSE) & b == 0, "Estimate"] = NA_real_
  DefOpt = options(contrasts=c("contr.SAS", "contr.SAS"))
  coef1 = coef1[colnames(model.matrix(Formula, Data)), , drop=FALSE]
  options(DefOpt)
  Res$coefficients = coef1
  Res$aliased = !is.numeric(coef1[,"Estimate"])

  Res$df = c(np, DFr, nc)
  Res$r.squared = 1 - SSE/SST

  if (DFr > 0) {
    Res$sigma = sqrt(MSE)
    Res$adj.r.squared = 1 - (1 - Res$r.squared) * (nr - fIntercept)/DFr
    Res$fstatistic = c(value=(SST - SSE)/(np - fIntercept)/MSE, numdf=(np - fIntercept), dendf=DFr)
  } else {
    Res$sigma = NaN
    Res$adj.r.squared = NaN
    Res$fstatistic = c(NaN, numdf=(np - fIntercept), dendf=DFr)
  }
  class(Res) = "summary.lm"
  return(Res)
}

