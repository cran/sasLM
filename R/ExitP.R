ExitP = function(Theta, bi, ti=NULL)
{
  K = length(bi)
  if (is.null(ti)) ti = (1:K)/K

  c0 = sqrt(outer(ti, ti, "/"))
  c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]

  ep = rep(NA, K) # cumulative exit probability
  ep[1] = (1 - pnorm(bi[1], mean=Theta*sqrt(ti[1]), sd=1))
  mi = Theta*sqrt(ti)

# pmvt (with noncentrality) is better than pmvnorm in calculating power and sample size. 
# But, Lan-DeMets used multi-variate normal rather than multi-variate noncentral t distributionh.
# I followed Lan-DeMets for the consistency with previous results.
  if (K < 11) {
    for (i in 2:K) ep[i] = 1 - pmvnorm(upper=bi[1:i], mean=mi[1:i], corr=c0[1:i, 1:i], algorithm=Miwa)
  } else {
    for (i in 2:K) ep[i] = 1 - pmvnorm(upper=bi[1:i], mean=mi[1:i], corr=c0[1:i, 1:i], seed=5) 
  }
  Res = cbind(exit.p = ep - c(0, ep[1:(K - 1)]), cum.exit.p = ep)
  attr(Res, "drift") = Theta
  return(Res)
}
