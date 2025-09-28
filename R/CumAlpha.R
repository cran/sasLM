CumAlpha = function(z, side=2, ti=NULL, c0=NULL, Seed=5)
{
  K = length(z)
  if (is.null(ti)) ti = (1:K)/K
  if (is.null(c0)) {
    c0 = sqrt(outer(ti, ti, "/"))
    c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]
  }
  pz = rep(NA, K)
  pz[1] = side*(1 - pnorm(z[1]))
  if (K < 11) { # Too slow with larger K (11 - 20)
    for (i in 2:K) pz[i] = side*(1 - pmvnorm(upper=z[1:i], corr=c0[1:i, 1:i], algorithm=Miwa))
  } else {
    for (i in 2:K) pz[i] = side*(1 - pmvnorm(upper=z[1:i], corr=c0[1:i, 1:i], seed=Seed))
  }
  return(cbind(ti, cum.alpha=pz))
}
