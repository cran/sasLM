ORmn = function(d0, conf.level=0.95, weight="MN", eps=1e-8)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(c(n1, n2) == 0)) stop("Check the input!")
  v0 = qchisq(conf.level, 1)
  p1 = y1/n1
  p2 = y2/n2
  w = n1*n2/(n1 + n2) # MH weight
  OR = sum(w/sum(w)*(p1/(1 - p1)/(p2/(1 - p2))))

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

  wor = function(or, eps=1e-8, MaxIter=50) {
    pw = n1*n2/(n1 + n2)
    for (i in 1:MaxIter) {
      p = p1p2or(or)
      r1 = p[, 1]
      r2 = p[, 2]
      w = n1*n2*((1 - r1)*r2)^2/(n1*r1*(1 - r1) + n2*r2*(1 - r2))
      if (sum(abs(w - pw)) < eps) break
      pw = w
    }
    return(w)
  }

  Obj = function(or) {
    p = p1p2or(or)
    r1 = p[, 1]
    r2 = p[, 2]
    if (weight == "MN") w = wor(or)    
    nume = sum(w*((p1 - r1)/(r1*(1 - r1)) - (p2 - r2)/(r2*(1 - r2))))
    v = (1/(n1*r1*(1 - r1)) + 1/(n2*r2*(1 - r2)))*(n1 + n2)/(n1 + n2 - 1)
#    return((nume^2/sum(w*w*v) - v0)^2)  # for nlminb
    return(nume^2/sum(w*w*v) - v0) # for uniroot
  }

  options(warn=-1)
#  lower = nlminb(OR, Obj, lower=0+eps, upper=OR-eps)$par
#  upper = nlminb(OR, Obj, lower=OR+eps)$par
  if (OR < eps) lower = 0
  else lower = uniroot(Obj, interval=c(eps, OR - eps))$root
  upper = uniroot(Obj, interval=c(OR + eps, 1e9))$root
  options(warn=1)

  nr = nrow(d0)
  Res1 = ORmn1(y1[1], n1[1], y2[1], n2[1], conf.level=conf.level)
  rownames(Res1) = NULL
  if (nr > 1) {
    for (i in 2:nr) {
      Res1 = rbind(Res1, ORmn1(y1[i], n1[i], y2[i], n2[i], conf.level=conf.level))
    }
    o1 = sum(w/sum(w)*p1/(1 - p1))
    o2 = sum(w/sum(w)*p2/(1 - p2))
    Res2 = c(odd1=o1, odd2=o2, OR=OR, lower=lower, upper=upper)
    Res1 = list(Strata=Res1, Common=Res2)
  }
  return(Res1)
}
