OBFBound = function(K, alpha=0.05, side=2, ti=NULL, c0=NULL) 
{
  if (is.null(ti)) ti = (1:K)/K
  iti = 1/ti
  if (is.null(c0)) {
    c0 = sqrt(outer(ti, ti, "/"))
    c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]
  }  
  cOBF = function(k) {
    zs = k*sqrt(iti)
    CumAlpha(zs, side=side)[K, "cum.alpha"] - alpha
  }
  cOBF = uniroot(cOBF, c(1, 10))$root
  z = cOBF*sqrt(iti)
  cum.alpha = CumAlpha(z, side=side)
  return(cbind(z, cum.alpha)[, c("ti", "z", "cum.alpha")])
}
