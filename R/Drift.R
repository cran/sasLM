Drift = function(bi, ti=NULL, Power=0.9)
{
  K = length(bi)
  if (is.null(ti)) ti = (1:K)/K
  c0 = sqrt(outer(ti, ti, "/"))
  c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]

# pmvt (with noncentrality) is better than pmvnorm in calculating power and sample size. 
# But, Lan-DeMets used multi-variate normal rather than multi-variate noncentral t distributionh.
# I followed Lan-DeMets for the consistency with previous results.
  if (K < 11) {
    fx = function(TH) 1 - pmvnorm(upper=bi, mean=TH*sqrt(ti), corr=c0, algorithm=Miwa) - Power
  } else {
    fx = function(TH) 1 - pmvnorm(upper=bi, mean=TH*sqrt(ti), corr=c0, seed=5) - Power
  }
  return(uniroot(fx, interval=c(0, 10))$root)
}
