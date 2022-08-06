RRmn1 = function(y1, n1, y2, n2, conf.level=0.95)
{
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(c(n1, n2) == 0)) stop("Check the input!")
  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  RR = p1/p2                 # point estimate of relative risk (RR)
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(rr) {  # find rr points of increased obj fx value (ofv) by v0
    mLL = function(p2d) -(log(dbinom(y1, n1, p2d*rr)) + log(dbinom(y2, n2, p2d)))
    p2d = nlminb(p2, mLL, lower=0, upper=1)$par
    p1d = p2d*rr                           # MLE p1d, p2d with fixed ratio rr
    var0 = (p1d*(1 - p1d)/n1 + RR*RR*p2d*(1 - p2d)/n2)*(n1 + n2)/(n1 + n2 - 1)
    return(((p1d - p2d*RR)^2/var0 - v0)^2) # find the root of increased ofv by v0
  }

  options(warn=-1)
  LL = nlminb(RR, Obj, lower=0, upper=RR)$par
  UL = nlminb(RR, Obj, lower=RR)$par
  options(warn=1)

  return(data.frame(p1 = p1, p2 = p2, RR = RR, lower = LL, upper = UL))
}
