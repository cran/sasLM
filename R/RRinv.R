RRinv = function(d0, conf.level=0.95)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")

  r1 = RR(y1, n1, y2, n2, conf.level=conf.level)
  if (any(r1$SElog < 1e-8)) warning("Note that standard error is too small!")
  
  thi = log(r1$RR)
  wi = 1/(r1$SElog)^2
  wi0 = y1*n1/(n1 + n2)
  pwi = wi0/sum(wi0)*100
  r1$pwi = pwi
  sumwi = sum(wi) ; sumwi

  th.hat = sum(wi*thi)/sumwi
  seth.hat = sqrt(1/sumwi)
  z.crit = qnorm(0.5 + conf.level/2)
  lower = exp(th.hat - z.crit*seth.hat) 
  upper = exp(th.hat + z.crit*seth.hat) 
  r2 = data.frame(RR = exp(th.hat), lower = lower, upper = upper)

  Q = sum(wi*(thi - th.hat)^2)
  k = length(y1) # number of strata
  pQ = 1 - pchisq(Q, k - 1)
  r3 = data.frame(Q=Q, prob=pQ)

  tau2 = (Q - (k - 1))/(sumwi - sum(wi^2)/sumwi) # method of moment
  wsi = 1/(1/wi + tau2)
  sumwsi = sum(wsi)
  pwsi = wsi/sumwsi*100
  r1$pwsi = pwsi

  th.hat.ran = sum(wsi*thi)/sumwsi
  se2 = sqrt(1/sumwsi)
  lower2 = exp(th.hat.ran - z.crit*se2)
  upper2 = exp(th.hat.ran + z.crit*se2)
  r4 = data.frame(RR = exp(th.hat.ran), lower = lower2, upper = upper2)
  Res = list(RRs = r1, Heterogeneity = r3, tau2 = tau2, Fixed = r2, Random = r4)
  return(Res)
}
