RDmn1 = function(y1, n1, y2, n2, conf.level=0.95, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  if (n1 > 1/eps | n2 > 1/eps) stop("Too large n1 or n2!")

  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  RD0 = p1 - p2               # point estimate of risk difference
  v0 = qchisq(conf.level, 1) # delta ofv for confidence interval

  Obj = function(rd) {  # find rd points of increased obj fx value (ofv) by v0
    L3 = n1 + n2                                   # eq 27, These could be float number!!!
    L2 = (n1 + 2*n2)*rd - L3 - y1 - y2
    L1 = (n2*rd - L3 - 2*y2)*rd + y1 + y2
    L0 = y2*rd*(1 - rd)
    q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3) # eq 28
    p = ifelse(abs(q) < eps, 0, sign(q)*sqrt(max(0, L2^2/(3*L3)^2 - L1/(3*L3))))
    a = (pi + ifelse(p == 0, acos(0), acos(min(1, max(-1, q/p^3)))))/3
    p2t = 2*p*cos(a) - L2/(3*L3)
    p1t = p2t + rd

    var0 = (p1t*(1 - p1t)/n1 + p2t*(1 - p2t)/n2)*L3/(L3 - 1)
    return((rd - RD0)^2/var0 - v0) # find roots of increased ofv by v0
  }

  options(warn=-1)
  if (RD0 < -1 + eps) {
    LL = -1
  } else {
    rTemp = try(uniroot(Obj, c(-1 + eps, RD0 - eps)), silent=T)
    if (!inherits(rTemp, "try-error")) {
      LL = rTemp$root
    } else {
      rTemp = try(uniroot(Obj, c(-1 + eps^2, RD0 - eps^2)), silent=T)
      if (!inherits(rTemp, "try-error")) { LL = rTemp$root
      } else { LL = NA }
    }
  }
  if (RD0 > 1 - eps) {
    UL = 1
  } else {
    rTemp = try(uniroot(Obj, c(RD0 + eps, 1 - eps)), silent=T)
    if (!inherits(rTemp, "try-error")) {
      UL = rTemp$root
    } else {
      rTemp = try(uniroot(Obj, c(RD0 + eps^2, 1 - eps^2)), silent=T)
      if (!inherits(rTemp, "try-error")) { UL = rTemp$root
      } else { UL = NA }
    }
  }
  options(warn=0)

  return(c(p1 = p1, p2 = p2, RD = RD0, lower = LL, upper = UL))
}
