PocockBound = function(K=2, alpha=0.05, side=2)
{
  fx = function(x) CumAlpha(x, K=K, side=side)[K, "cum.alpha"] - alpha
  Z0 = uniroot(fx, interval=c(0, 8))$root
  t1 = CumAlpha(Z0, K=K, side=side)
  attr(Z0, "Cumulative Alpha") = t1
  return(Z0)
}

