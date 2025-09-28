PocockBound = function(K=2, alpha=0.05, side=2)
{
  fx = function(x) CumAlpha(rep(x, K), side=side)[K, "cum.alpha"] - alpha
  z0 = uniroot(fx, interval=c(0, 8))$root
  t1 = CumAlpha(rep(z0, K), side=side)
  attr(z0, "Cumulative Alpha") = t1
  return(z0)
}
