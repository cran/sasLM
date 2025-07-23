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

  p1p2or = function(or) {
    if (or == 1) {
      p1t = (y1 + y2)/(n1 + n2)
      p2t = p1t
    } else {
      A = n2*(or - 1)                         # eq 13
      B = n1*or + n2 - (y1 + y2)*(or - 1)
      C1 = -(y1 + y2)
      p2t = (-B + sqrt(B*B - 4*A*C1))/(2*A)
      p1t = p2t*or/(1 + p2t*(or - 1))
    }
    return(cbind(p1t, p2t))
  }

  wor = function(or, MaxIter=50) {
    pw = n1*n2/(n1 + n2)
    for (i in 1:MaxIter) {
      p = p1p2or(or)
      r1 = p[, 1]
      r2 = p[, 2]
      w = n1*n2*((1 - r1)*r2)^2/(n1*r1*(1 - r1) + n2*r2*(1 - r2))
      if (sum(abs(w - pw)) < 1e-8) break
      pw = w
    }
    return(w)
  }

  w = n1*n2/(n1 + n2)  
  Obj = function(or) {
    p = p1p2or(or)
    r1 = p[, 1]
    r2 = p[, 2]
    w = wor(or)    
    w <<- w/sum(w)
    nume = sum(w*((p1 - r1)/(r1*(1 - r1)) - (p2 - r2)/(r2*(1 - r2))))
    v = (1/(n1*r1*(1 - r1)) + 1/(n2*r2*(1 - r2)))*(n1 + n2)/(n1 + n2 - 1)
    return(sqrt(nume^2/sum(w*w*v)) - z0)
  }

  OR0 = ORcmh(d0)$Common$OR[[1]]
  options(warn=-1)
  if (OR0 < eps) { 
    lower = 0
  } else { 
    lower = uniroot(Obj, interval=c(eps, OR0 - eps))$root 
  }
  upper = uniroot(Obj, interval=c(OR0 + eps, 1e9))$root
  options(warn=0)

  for (i in 2:nr) {
    Res1 = rbind(Res1, ORmn1(y1[i], n1[i], y2[i], n2[i], conf.level=conf.level))
  }
  rownames(Res1) = rownames(d0)    
  p1s = sum(w/sum(w)*p1)
  p2s = sum(w/sum(w)*p2)  
  o1 = p1s/(1 - p1s)
  o2 = p2s/(1 - p2s)
  Res2 = c(odd1=o1, odd2=o2, OR=o1/o2, lower=lower, upper=upper)
  Res1 = list(Strata=Res1, Common=Res2)

  return(Res1)
}
