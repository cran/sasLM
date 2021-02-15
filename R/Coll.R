Coll = function(Formula, Data)
{
  rx = lm(Formula, Data)
  np = length(rx$coefficients)

  d1 = model.frame(Formula, Data)[,-1]
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
  p1 = SVDz$v %*% diag(1/SVDz$d)
  p2 = t(p1^2)
  Res2 = prop.table(p2 %*% diag(rowSums(p2)), 2)

  ZpZ = crossprod(Z)
  ev = eigen(ZpZ/diag(ZpZ))$values
  cn = sqrt(ev[1]/ev)

  Res2 = cbind(ev, cn, Res2)
  colnames(Res2) = c("Eigenvalue", "Cond. Index", "(Intercept)", Names)

  Res = list(Res1, Res2)
  names(Res) = c("Tolerance and VIF", "Collinearity Diagnosis")

  return(Res)
}
