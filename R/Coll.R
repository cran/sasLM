Coll = function(Formula, Data)
{
  Terms = terms(Formula, data=Data)
  if (!attr(Terms, "response")) stop("Dependent variable should be provided!")
  rx = lm(Formula, Data)
  np = length(rx$coefficients)

  mf = model.frame(Formula, Data)
  for (i in 1:ncol(mf)) {
    if (!is.numeric(mf[, i])) stop("All variables should be numeric!")
  }

  d1 = mf[, -1, drop=FALSE]
  Names = colnames(d1)
  nName = NCOL(d1)
  Tol = rep(NA, nName)
  for (i in 1:nName) {
    f1 = as.formula(paste0('`', Names[i], '` ', "~ ."))
    r1 = lm(f1, d1)
    Tol[i] = 1 - (summary(r1)$r.squared)
  }
  VIF = 1/Tol

  Res1 = cbind(Tol, VIF)
  rownames(Res1) = Names

  Z = scale(model.matrix(Formula, Data), center=FALSE)
  SVDz = svd(Z)
  if (length(SVDz$d) == 1) {
    p1 = SVDz$v %*% diag(1/as.matrix(SVDz$d))
  } else {
    p1 = SVDz$v %*% diag(1/SVDz$d)
  }

  p2 = t(p1^2)
  if (NCOL(p2) == 1) {
    Res2 = prop.table(p2 %*% diag(as.matrix(rowSums(p2))), 2)
  } else {
    Res2 = prop.table(p2 %*% diag(rowSums(p2)), 2)
  }

  ZpZ = crossprod(Z)
  ev = eigen(ZpZ/diag(ZpZ))$values
  cn = sqrt(ev[1]/ev)

  Res2 = cbind(ev, cn, Res2)
  if (attr(Terms, "intercept")) {
    colnames(Res2) = c("Eigenvalue", "Cond. Index", "(Intercept)", Names)
  } else {
    colnames(Res2) = c("Eigenvalue", "Cond. Index", Names)    
  }

  Res = list(Res1, Res2)
  names(Res) = c("Tolerance and VIF", "Collinearity Diagnosis")

  return(Res)
}
