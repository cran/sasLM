Coll = function(Formula, Data)
{
  Terms = terms(Formula, data=Data)
  if (!attr(Terms, "response")) stop("Dependent variable should be provided!")
  rx = lm(Formula, Data)

  mf = model.frame(Formula, Data)
  for (i in 1:ncol(mf)) {
    if (!is.numeric(mf[, i])) stop("All variables should be numeric!")
  }

  X = model.matrix(rx)
  V = vcov(rx)
  if (colnames(X)[1] == "(Intercept)") {
    X = X[, -1, drop=F]
    V = V[-1, -1]
  }
  colNames = colnames(X)
  VIF = diag(solve(cov2cor(V)))
  Tol = 1/VIF

  Res1 = cbind(Tol, VIF)
  rownames(Res1) = colNames

  e = eigen(cor(X))
  ev= e$values
  names(ev) = colNames
  cn = sqrt(ev[1]/ev)
  pr = prop.table(t(e$vectors^2)/ev, margin=2)
  colnames(pr) = colNames
  
  Res2 = cbind(ev, cn, pr)
  colnames(Res2) = c("Eigenvalue", "Cond. Index", colNames)

  Res = list(Res1, Res2)
  names(Res) = c("Tolerance and VIF", "Collinearity Diagnostics")

  return(Res)
}
