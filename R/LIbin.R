LIbin = function(y, n, k, conf.level=0.95, eps=1e-8)
{
  if (length(y) != 1 | length(n) != 1) stop("Current version supports only length 1 data")
  if (any(y < 0) | any(n < 0) | any(!is.finite(y)) | any(!is.finite(n))) stop("Check the input!")
  logk = ifelse(missing(k), qf(conf.level, 1, max(n - 1, 1))/2, log(k))
  logk = min(logk, log(2/(1 - conf.level)))
  p0 = y/n
  if (y == 0 | y == n) {
    LE = logk/sqrt(n)                         # initial guess of the limit of error
    maxLogL = 0                               # lchoose removed
  } else {
    LE = logk*sqrt(p0*(1 - p0))               # initial guess of the limit of error
    maxLogL = y*log(p0) + (n - y)*log(1 - p0) # lchoose removed
  }

  Obj = function(p)  {
    logL = ifelse(p < 1e-300 | p > 1 - 1e-300, 0, y*log(p) + (n - y)*log(1 - p))
    return(maxLogL - logL - logk)
  }

  if (p0 < eps) {
    LL = 0
  } else {
    rTemp = try(uniroot(Obj, c(max(eps, p0 - LE), p0 + eps)), silent=TRUE)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = max(p0 - LE, 0) }
  }

  if (p0 > 1 - eps) {
    UL = 1
  } else {
    rTemp = try(uniroot(Obj, c(p0 - eps, min(p0 + LE, 1 - eps))), silent=TRUE)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = min(p0 + LE, 1) }
  }

  Res = c(y = y, n = n, PE=p0, LL=LL, UL=UL, k=exp(logk))
  return(Res)
}
