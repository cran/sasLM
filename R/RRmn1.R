RRmn1 = function(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  RR0 = p1/p2                # point estimate of relative risk (RR)
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj1 = function(rr) {
    A = (n1 + n2)*rr
    B = n1*rr + y1 + n2 + y2*rr
    C1 = y1 + y2
    p1t = (B - sqrt(B*B - 4*A*C1))/(2*A)
    p2t = p1t*rr

    var0 = n1*n2*p2t/(n1*(rr - p2t) + n2*(1 - p2t))
    chisq0 = ((y1 - n1*p2t)/(1 - p2t))^2/var0
    return(chisq0 - v0)
  }

  Obj2 = function(rr) { # Test and control group switched
    A = (n1 + n2)*rr
    B = n2*rr + y2 + n1 + y1*rr
    C1 = y1 + y2
    p1t = (B - sqrt(B*B - 4*A*C1))/(2*A)
    p2t = p1t*rr

    var0 = n1*n2*p2t/(n2*(rr - p2t) + n1*(1 - p2t))
    chisq0 = ((y2 - n2*p2t)/(1 - p2t))^2/var0
    return(chisq0 - v0)
  }

  if (y1 == 0 & y2 == 0) {                         # Case 1: p1 = 0 & p2 = 0
    LL = 0
    UL = Inf
  } else if (y1 == 0 & y2 > 0 & y2 < n2) {         # Case 2: p1 = 0 & 0 < p2 < 1
    LL = 0
    UL = uniroot(Obj1, c(eps, 1e9))$root
  } else if (y1 == 0 & y2 == n2) {                 # Case 3: p1 = 0 & p2 = 1
    LL = 0
    UL = uniroot(Obj1, c(eps, 1e9))$root
  } else if (y1 > 0 & y1 < n1 & y2 == 0) {         # Case 4: 0 < p1 < 1 & p2 = 0
    LL = uniroot(Obj1, c(eps, 1e9))$root
    UL = Inf
  } else if (y1 > 0 & y1 < n1 & y2 > 0 & y2 < n2) {# Case 5: 0 < p1 < 1 & 0 < p2 < 1
    LL = uniroot(Obj1, c(eps, RR0))$root
    UL = uniroot(Obj1, c(RR0, 1e9))$root
  } else if (y1 > 0 & y1 < n1 & y2 == n2) {        # Case 6: 0 < p1 < 1 & p2 = 1
    LL = uniroot(Obj1, c(eps, RR0))$root
    UL = uniroot(Obj1, c(RR0, 1e9))$root
  } else if (y1 == n1 & y2 == 0) {                 # Case 7: p1 = 1 &  p2 = 0
    LL = 1/uniroot(Obj2, c(eps, 1e9))$root
    UL = Inf
  } else if (y1 == n1 & y2 > 0 & y2 < n2) {        # Case 8: p1 = 1 & 0 < p2 < 1
    LL = 1/uniroot(Obj2, c(1/RR0, 1e9))$root
    UL = 1/uniroot(Obj2, c(eps, 1/RR0))$root
  } else if (y1 == n1 & y2 == n2) {                # Case 9: p1 = 1 &  p2 = 1
    LL = n1/(n1 + v0)
    UL = (n2 + v0)/n2
  }

  return(c(p1 = p1, p2 = p2, RR = RR0, lower = LL, upper = UL))
}
