RDmn1 = function(y1, n1, y2, n2, conf.level=0.95)
{
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(c(n1, n2) == 0)) stop("Check the input!")
  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  dp0 = p1 - p2              # point estimate of risk difference (RD)
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(dp) {  # find dp points of increased obj fx value (ofv) by v0
    mLL = function(p1d) -(log(dbinom(y1, n1, p1d)) + log(dbinom(y2, n2, p1d - dp)))
    p1d = nlminb(p1, mLL, lower=max(0, dp0), upper=min(1, dp0 + 1))$par
    p2d = max(0, min(p1d - dp, 1)) # MLE p1d, p2d with fixed delta p (dp)
    var0 = (p1d*(1 - p1d)/n1 + p2d*(1 - p2d)/n2)*(n1 + n2)/(n1 + n2 - 1)
    return(((dp - dp0)^2/var0 - v0)^2) # find roots of increased ofv by v0
  }

  options(warn=-1)
  LL = nlminb(dp0, Obj, lower=max(-1, dp0 - 1), upper=dp0)$par
  UL = nlminb(dp0, Obj, lower=dp0, upper=min(dp0 + 1, 1))$par
  options(warn=1)

  return(data.frame(p1 = p1, p2 = p2, RD = dp0, lower = LL, upper = UL))
}
