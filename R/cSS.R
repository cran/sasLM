cSS = function(K, rx, eps=1e-8)
{
  b = rx$coefficients
  iiv = G2SWEEP(t(K) %*% rx$g2 %*% K, Augmented=FALSE, eps=eps)
  Q = t(t(K) %*% b) %*% iiv %*% t(K) %*% b
  Df = attr(iiv, "rank")
  MS = ifelse(Df > 0, Q/Df, NA)
  if (rx$DFr > 0) {
    Fval = MS/(rx$SSE/rx$DFr)
    Pval = 1 - pf(Fval, Df, rx$DFr)
  } else {
    Fval = NA
    Pval = NA
  }
  Result = c(Df=Df, 'Sum Sq'=Q, 'Mean Sq'=MS, 'F value'=Fval, 'Pr(>F)'=Pval)
  return(Result)
}
