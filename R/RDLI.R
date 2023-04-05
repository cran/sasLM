RDLI = function(y1, n1, y2, n2, conf.level=0.95, k, eps=1e-8)
{
  if (length(y1) != 1 | length(n1) != 1 | (y1 < 0) | (n1 < 0) | !is.finite(y1) | !is.finite(n1)) stop("Check the input!")
  if (length(y2) != 1 | length(n2) != 1 | (y2 < 0) | (n2 < 0) | !is.finite(y2) | !is.finite(n2)) stop("Check the input!")
  if (n1 > 1/eps | n2 > 1/eps) stop("Too large n1 or n2!")

  p1 = y1/n1
  p2 = y2/n2
  RD0 = p1 - p2
  if (p1 == 0) p1 = eps # to prevent likelihood = 0 (loglihood is -inf)
  if (p2 == 0) p2 = eps # to prevent likelihood = 0 (loglihood is -inf)
  if (p1 == 1) p1 = 1 - eps # to prevent likelihood = 0 (loglihood is -inf)
  if (p2 == 1) p2 = 1 - eps # to prevent likelihood = 0 (loglihood is -inf)

  maxLL = y1*log(p1) + (n1 - y1)*log(1 - p1) + y2*log(p2) + (n2 - y2)*log(1 - p2) # lchoose part removed
  maxLL = ifelse(is.finite(maxLL), maxLL, 0)

  logk = ifelse(missing(k), qf(conf.level, 1, n1 + n2 - 1)/2, log(k))
  logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05

  Obj = function(rd) {
    L3 = n1 + n2                                   # eq 27, These could be float number!!!
    L2 = (n1 + 2*n2)*rd - L3 - y1 - y2
    L1 = (n2*rd - L3 - 2*y2)*rd + y1 + y2
    L0 = y2*rd*(1 - rd)
    q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3) # eq 28
    if (!is.finite(q) | abs(q) < eps) q = 0
    p = sign(q)*sqrt(max(0, L2^2/(3*L3)^2 - L1/(3*L3)))
    a = (pi + ifelse(p == 0, acos(0), acos(min(max(q/p^3, -1), 1))))/3
    p2t = min(max(eps, 2*p*cos(a) - L2/(3*L3)), 1 - eps)
    p1t = min(max(eps, p2t + rd), 1 - eps)

    ll = y1*log(p1t) + (n1 - y1)*log(1 - p1t) + y2*log(p2t) + (n2 - y2)*log(1 - p2t) # lchoose part removed
    ll = ifelse(is.finite(ll), ll, 0)
    return(maxLL - ll - logk)
  }
  options(warn=-1)
  if (RD0 < -1 + eps) {
    LL = -1
  } else {
    rTemp = try(uniroot(Obj, c(-1 + eps, RD0 + eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = -1 }
  }
  if (RD0 > 1 - eps) {
    UL = 1
  } else {
    rTemp = try(uniroot(Obj, c(RD0 - eps, 1 - eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = 1 }
 }
  options(warn=0)

  Res = c(p1 = y1/n1, p2 = y2/n2, RD = RD0, lower = LL, upper = UL, k = exp(logk))
  return(Res)
}
