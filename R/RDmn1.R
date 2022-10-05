RDmn1 = function(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  RD = p1 - p2               # point estimate of risk difference
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(rd) {  # find rd points of increased obj fx value (ofv) by v0
#    mLL = function(p1t) -dbinom(y1, n1, p1t, log=T) + dbinom(y2, n2, p1t - rd, log=T)
#    p1t = nlminb(p1, mLL, lower=max(0, RD), upper=min(1, RD + 1))$par
#    p2t = max(0, min(p1t - rd, 1)) # MLE p1t, p2t with fixed delta p (rd)
    L3 = n1 + n2                                   # eq 27
    L2 = (n1 + 2*n2)*rd - L3 - y1 - y2
    L1 = (n2*rd - L3 - 2*y2)*rd + y1 + y2
    L0 = y2*rd*(1 - rd)
    q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3) # eq 28
    p = sign(q)*sqrt(L2^2/(3*L3)^2 - L1/(3*L3))
    a = (pi + acos(q/p^3))/3
    p2t = 2*p*cos(a) - L2/(3*L3)
    p1t = p2t + rd

    var0 = (p1t*(1 - p1t)/n1 + p2t*(1 - p2t)/n2)*(n1 + n2)/(n1 + n2 - 1)
    return(((rd - RD)^2/var0 - v0)^2) # find roots of increased ofv by v0
  }

  options(warn=-1)
  LL = nlminb(max(-1 + eps, RD - eps), Obj, lower=max(-1, RD - 1), upper=RD)$par
  UL = nlminb(min(1 - eps, RD + eps), Obj, lower=RD, upper=min(RD + 1, 1))$par
  options(warn=1)

  return(c(p1 = p1, p2 = p2, RD = RD, lower = LL, upper = UL))
}
