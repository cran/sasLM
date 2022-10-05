RDinv = function(d0, conf.level=0.95)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")

  r1 = RD(y1, n1, y2, n2, conf.level=conf.level)
  wi = 1/(r1$SE)^2
  sumwi = sum(wi)
  pe = sum(wi*r1$RD)/sumwi
  se = sqrt(1/sumwi)
  z.crit = qnorm(0.5 + conf.level/2)
  lower = pe - z.crit*se
  upper = pe + z.crit*se
  r2 = data.frame(PE=pe, SE=se, lower=lower, upper=upper)
  Q = sum(wi*(r1$RD - pe)^2)
  
  k = length(wi)
  pQ = 1 - pchisq(Q, k - 1)
  r3 = data.frame(Q=Q, prob=pQ)
  
  tau2 = max(0, (Q - (k - 1))/(sumwi - sum(wi^2)/sumwi))
  wsi = 1/(1/wi + tau2)
  sumwsi = sum(wsi)
  pe2 = sum(wsi*r1$RD)/sumwsi
  se2 = sqrt(1/sumwsi)
  lower2 = pe2 - z.crit*se2
  upper2 = pe2 + z.crit*se2
  r4 = data.frame(PE=pe2, SE=se2, lower=lower2, upper=upper2)
  
  Res = list(RDs = r1, Heterogeneity = r3, tau2 = tau2, Fixed = r2, Random = r4)
  return(Res)
}
