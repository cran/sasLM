regD = function(formula, data)
{
  X = model.matrix(formula, data)
  y = model.frame(formula, data)[,1]

  n = nrow(X)
  p = qr(X)$rank
  if (p != ncol(X)) stop("Model matrix is not full rank. Consider other way!")

  VIF = solve(crossprod(X))        # Variance Inflation Factor
  b = VIF %*% crossprod(X, y)      # Beta Hat
  yhat = X %*% b                   # Predicted
  e = y - yhat                     # Residual
  SSE = sum(e*e)                   # e'e or as.numeric(crossprod(e))
  DFr = n - p                      # Degree of freedom
  MSE = SSE/DFr
  bVar = VIF * MSE                 # Variance-Covariance Matrix of Beta Hat
  bSE = sqrt(diag(bVar))           # Standard Error of Beta Hat

  h = diag(X %*% VIF %*% t(X))     # hat
  Rse = sqrt(MSE*(1 - h))          # SE of RStudent residual
  Rst = e/Rse                      # RStudent residual
  CooksD = e*e/(p*MSE)*h/(1 - h)^2 # CooksD

  MSEi = (SSE - e*e/(1 - h))/(DFr - 1)
  sdr = e/sqrt(MSEi*(1 - h))       # Studentized Deleted Residual
  DFFITS = sdr*sqrt(h/(1 - h))

  COVRATIO = matrix(nrow=n)
  dCov = det(bVar)
  bi = matrix(nrow=n, ncol=p)
  for (i in 1:n) {
    Xmi = X[-i, ]
    iXpXmi = solve(crossprod(Xmi))
    COVRATIO[i] = det(MSEi[i] * iXpXmi)/dCov
    bi[i, ] = iXpXmi %*% crossprod(Xmi, y[-i])
  }
  bm = matrix(rep(t(b), n), ncol=p, byrow=TRUE)
  DFBETAS = (bm - bi)/sqrt(MSEi %*% diag(VIF))
  colnames(DFBETAS) = rownames(b)

  tVal = b/bSE
  pVal = 2*pt(-abs(tVal), DFr)
  Res0 = cbind(b, bSE, tVal, pVal)
  colnames(Res0) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  class(Res0) = "anova"

  Res1 = cbind(yhat, e, Rse, Rst, h, CooksD, sdr, DFFITS, COVRATIO)
  colnames(Res1) = c("Predicted", "Residual", "se_resid", "RStudent", "Leverage", "Cook's D", "sdResid", "DFFITS", "COVRATIO")

  Res = list(Res0, Res1, DFBETAS)
  names(Res) = c("Coefficients", "Diagnostics", "DFBETAS")
  return(Res)
}
