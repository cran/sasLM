ORmn = function(d0, conf.level=0.95, eps=1e-8)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")
  z0 = qnorm(0.5 + conf.level/2)
  p1 = y1/n1
  p2 = y2/n2

  nr = nrow(d0)
  Res1 = ORmn1(y1[1], n1[1], y2[1], n2[1], conf.level=conf.level)
  rownames(Res1) = NULL
  if (nr == 1) return(Res1)

  lambda = (n1 + n2)/(n1 + n2 - 1) # variance bias correction as Miettinen-Nurminen

  p1p2or = function(or) { # restricted MLEs under common odds ratio
    if (or == 1) {
      p2t = (y1 + y2)/(n1 + n2)
      p1t = p2t
    } else {
      A = n2*(or - 1)                         # eq 13
      B = n1*or + n2 - (y1 + y2)*(or - 1)
      C1 = -(y1 + y2)
      p2t = (-B + sqrt(B*B - 4*A*C1))/(2*A)
      p1t = p2t*or/(1 + p2t*(or - 1))
    }
    p1t = pmin(pmax(p1t, eps), 1 - eps)
    p2t = pmin(pmax(p2t, eps), 1 - eps)
    return(cbind(p1t, p2t))
  }

  wor = function(or) { # inverse variance weight at restricted MLEs
    p = p1p2or(or)
    r1 = p[, 1]
    r2 = p[, 2]
    return(1/(1/(n1*r1*(1 - r1)) + 1/(n2*r2*(1 - r2))))
  }

  zScore = function(or) { # inverse variance weighted score statistic
    p = p1p2or(or)
    r1 = p[, 1]
    r2 = p[, 2]
    s = (p1 - r1)/(r1*(1 - r1)) - (p2 - r2)/(r2*(1 - r2))
    s = s - (r1 - r2)/(n1*r1*(1 - r1) + n2*r2*(1 - r2)) # score bias correction as Laud
    s[(y1 == 0 & y2 == 0) | (y1 == n1 & y2 == n2)] = 0
    w = 1/(1/(n1*r1*(1 - r1)) + 1/(n2*r2*(1 - r2)))
    return(sum(w/sum(w)*s)/sqrt(sum(lambda*w)/sum(w)^2))
  }

  OR0 = suppressWarnings(uniroot(zScore, interval=c(eps, 1e9), tol=1e-10)$root)
  if (zScore(eps) - z0 < 0) {
    lower = 0
  } else {
    lower = suppressWarnings(uniroot(function(or) zScore(or) - z0, interval=c(eps, OR0), tol=1e-10)$root)
  }
  upper = suppressWarnings(uniroot(function(or) zScore(or) + z0, interval=c(OR0, 1e9), tol=1e-10)$root)

  for (i in 2:nr) {
    Res1 = rbind(Res1, ORmn1(y1[i], n1[i], y2[i], n2[i], conf.level=conf.level))
  }
  rownames(Res1) = rownames(d0)
  w = wor(OR0)
  p1s = sum(w/sum(w)*p1)
  p2s = sum(w/sum(w)*p2)
  o1 = p1s/(1 - p1s)
  o2 = p2s/(1 - p2s)
  Res2 = c(odd1=o1, odd2=o2, OR=OR0, lower=lower, upper=upper)
  Res1 = list(Strata=Res1, Common=Res2)

  return(Res1)
}
