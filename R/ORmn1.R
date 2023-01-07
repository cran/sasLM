ORmn1 = function(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  if (any((n1 - y1)*(n2 - y2) == 0)) stop("OR is 0 or infinity or NaN(not a number)! Consider substracting a tiny number from the numerator.")
  p1 = y1/n1
  p2 = y2/n2
  o1 = y1/(n1 - y1)          # odd of test (active) group
  o2 = y2/(n2 - y2)          # odd of control (placebo) group
  OR = o1/o2                 # point estimate of odds ratio (OR)
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(or) {  # find or points of increased obj fx value (ofv) by v0
#    mLL = function(p2d) {
#      p1d = p2d*or/(1 + p2d*(or - 1))
#      -dbinom(y1, n1, p1d, log=T) + dbinom(y2, n2, p2d, log=T)
#    }
#    p2d = nlminb(p2, mLL, lower=0, upper=1)$par
#    p1d = p2d*or/(1 + p2d*(or - 1))   # MLE p1d, p2d with fixed odds ratio or
    A = n2*(or - 1)                         # eq 13
    B = n1*or + n2 - (y1 + y2)*(or - 1)
    C1 = -(y1 + y2)
    p2t = (-B + sqrt(B*B - 4*A*C1))/(2*A)
    p1t = p2t*or/(1 + p2t*(or - 1))

    nume = (p1 - p1t)/(p1t*(1 - p1t)) - (p2 - p2t)/(p2t*(1 - p2t))
    var0 = (1/(n1*p1t*(1 - p1t)) + 1/(n2*p2t*(1 - p2t)))*(n1 + n2)/(n1 + n2 - 1)
    return((nume^2/var0 - v0)^2)       # find the root of increased ofv by v0
  }

  options(warn=-1)
  LL = nlminb(max(eps, OR - eps), Obj, lower=0, upper=OR)$par
  UL = nlminb(OR + eps, Obj, lower=OR)$par
  options(warn=1)

  return(c(odd1 = o1, odd2 = o2, OR = OR, lower = LL, upper = UL))
}
