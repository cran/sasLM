CumAlpha = function(x, K=2, side=2)
{
  ti = (1:K)/K
  c0 = sqrt(outer(ti, ti, "/"))
  c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]

  px = rep(NA, K)
  px[1] = side*(1 - pnorm(x))
  if (K < 11) { # Too slow with larger K (11 - 20)
    for (i in 2:K) px[i] = side*(1 - pmvnorm(upper=rep(x, i), corr=c0[1:i, 1:i], algorithm=Miwa))
  } else {
    for (i in 2:K) px[i] = side*(1 - pmvnorm(upper=rep(x, i), corr=c0[1:i, 1:i], seed=5))
  }

  return(cbind(ti, cum.alpha=px))
}
