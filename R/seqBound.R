seqBound = function(ti, alpha=0.05, side=2, t2=NULL, asf=1)
{
  K = length(ti)
  alpha0 = alpha/side
  if (asf == 1) {        # O'Brien-Flemming
    ai = 2*(1 - pnorm(qnorm(1 - alpha0/2)/sqrt(ti)))
  } else if (asf == 2) { # Pocock
    ai = alpha0*log(1 + 1.718281828459*ti)
  } else if (asf == 3) {
    ai = alpha0*ti
  } else if (asf == 4) {
    ai = alpha0*ti^1.5
  } else if (asf == 5) {
    ai = alpha0*ti^2
  } else {
    stop("Unknown alpha spending function")
  }

  if (is.null(t2)) t2 = ti
  c0 = sqrt(outer(t2, t2, "/"))
  c0[lower.tri(c0)] = t(c0)[lower.tri(c0)]

  bi = rep(NA, K)
  bi[1] = qnorm(1 - ai[1])
  px = rep(NA, K)
  px[1] = side*(1 - pnorm(bi[1]))

  if (K < 11) { # too slow with lareger K (12 - 20)
    fx = function(qx) 1 - pmvnorm(upper=c(bi[1:(i - 1)], qx), corr=c0[1:i, 1:i], algorithm=Miwa) - ai[i]
    for (i in 2:K) bi[i] = uniroot(fx, interval=c(0, 8))$root
    for (i in 2:K) px[i] = side*(1 - pmvnorm(upper=bi[1:i], corr=c0[1:i, 1:i], algorithm=Miwa))
  } else {
    fx = function(qx) 1 - pmvnorm(upper=c(bi[1:(i - 1)], qx), corr=c0[1:i, 1:i], seed=5) - ai[i]
    for (i in 2:K) bi[i] = uniroot(fx, interval=c(0, 8))$root
    for (i in 2:K) px[i] = side*(1 - pmvnorm(upper=bi[1:i], corr=c0[1:i, 1:i], seed=5))
  }

  return(cbind(time=ti, up.bound=bi, cum.alpha=px))
}
