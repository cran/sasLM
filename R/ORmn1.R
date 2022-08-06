ORmn1 = function(y1, n1, y2, n2, conf.level=0.95)
{
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(c(n1, n2) == 0)) stop("Check the input!")
  if (any(c(n1 - y1, n2 - y2) == 0)) stop("Check the input!")
  p1 = y1/n1
  p2 = y2/n2
  o1 = y1/(n1 - y1)          # odd of test (active) group
  o2 = y2/(n2 - y2)          # odd of control (placebo) group
  OR = o1/o2                 # point estimate of odds ratio (OR)
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(or) {  # find or points of increased obj fx value (ofv) by v0
    mLL = function(p2d) {
      p1d = p2d*or/(1 + p2d*(or - 1))
      -(log(dbinom(y1, n1, p1d)) + log(dbinom(y2, n2, p2d)))
    }
    p2d = nlminb(p2, mLL, lower=0, upper=1)$par
    p1d = p2d*or/(1 + p2d*(or - 1))   # MLE p1d, p2d with fixed odds ratio or
    div = (p1 - p1d)/(p1d*(1 - p1d)) - (p2 - p2d)/(p2d*(1 - p2d))
    var0 = (1/(n1*p1d*(1 - p1d)) + 1/(n2*p2d*(1 - p2d)))*(n1 + n2)/(n1 + n2 - 1)
    return((div^2/var0 - v0)^2)       # find the root of increased ofv by v0
  }

  options(warn=-1)
  LL = nlminb(OR, Obj, lower=0, upper=OR)$par
  UL = nlminb(OR, Obj, lower=OR)$par
  options(warn=1)

  return(data.frame(odd1 = o1, odd2 = o2, OR = OR, lower = LL, upper = UL))
}
