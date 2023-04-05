ORmn1 = function(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  p1 = y1/n1
  p2 = y2/n2
  o1 = y1/(n1 - y1)                 # odd of test (active) group
  o2 = y2/(n2 - y2)                 # odd of control (placebo) group
  OR0 = y1/y2*(n2 - y2)/(n1 - y1)   # point estimate of odds ratio (OR)
  v0 = qchisq(conf.level, 1)        # delta ofv for confidence interval

  Obj = function(or) {  # find or points of increased obj fx value (ofv) by v0
    A = n2*(or - 1)                         # eq 13
    B = n1*or + n2 - (y1 + y2)*(or - 1)
    C1 = -(y1 + y2)
    p2t = (-B + sqrt(B*B - 4*A*C1))/(2*A)
    p1t = p2t*or/(1 + p2t*(or - 1))

    nume = (p1 - p1t)/(p1t*(1 - p1t)) - (p2 - p2t)/(p2t*(1 - p2t))
    var0 = (1/(n1*p1t*(1 - p1t)) + 1/(n2*p2t*(1 - p2t)))*(n1 + n2)/(n1 + n2 - 1)
    return(nume^2/var0 - v0)         # find the root of increased ofv by v0
  }

  if (y1 == 0 & y2 == 0) {           # Case 1; o1=0 & o2=0
    LL = 0
    UL = Inf
  } else if (y1 == 0 & y2 < n2) {    # Case 2: o1=0 & 0 < o2 < Inf
    LL = 0
    UL = uniroot(Obj, c(eps, 1e9))$root
  } else if (y1 == 0 & y2 == n2) {   # Case 3: o1=0 & o2=Inf 
    LL = 0
    UL = uniroot(Obj, c(eps, 1e9))$root
  } else if (y1 < n1 & y2 == 0) {    # Case 4: o1=Inf & o2=0
    LL = uniroot(Obj, c(eps, 1e9))$root
    UL = Inf
  } else if (y1 < n1 & y2 < n2) {    # Case 5: 0 < o1 < Inf & 0 < o2 < Inf
    LL = uniroot(Obj, c(eps, OR0 + eps))$root
    UL = uniroot(Obj, c(OR0 - eps, 1e9))$root
  } else if (y1 < n1 & y2 == n2) {   # Case 6: 0 < o1 < 1 & o2=Inf
    LL = 0
    UL = uniroot(Obj, c(eps, 1e9))$root
  } else if (y1 == n1 & y2 == 0) {   # Case 7: o1=Inf & o2=0
    LL = uniroot(Obj, c(eps, 1e9))$root
    UL = Inf
  } else if (y1 == n1 & y2 < n2) {   # Case 8: o1=Inf & 0 < o2 < Inf
    LL = uniroot(Obj, c(eps, 1e9))$root
    UL = Inf
  } else if (y1 == n1 & y2 == n2) {  # Case 9: o1=Inf & o2=Inf
    LL = 0
    UL = Inf
  }

  return(c(odd1 = o1, odd2 = o2, OR = OR0, lower = LL, upper = UL))
}

