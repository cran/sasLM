RDmn = function(d0, conf.level=0.95, eps=1e-8)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]

  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")
  if (any(n1 > 1/eps) | any(n2 > 1/eps)) stop("Too large n1 or n2!")

  Res1 = RDmn1(y1[1], n1[1], y2[1], n2[1], conf.level=conf.level)
  rownames(Res1) = NULL

  nr = nrow(d0)
  if (nr == 1) return(Res1)
  
  v0 = qchisq(conf.level, 1)
  p1 = y1/n1
  p2 = y2/n2

  p1p2rd = function(rd) {
    L3 = n1 + n2
    L2 = (n1 + 2*n2)*rd - L3 - y1 - y2
    L1 = (n2*rd - L3 - 2*y2)*rd + y1 + y2
    L0 = y2*rd*(1 - rd)
    q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3) # eq 28
    p = sign(q)*sqrt(max(0, L2^2/(3*L3)^2 - L1/(3*L3)))
    a = (pi + ifelse(p == 0, acos(0), acos(min(1, max(-1, q/p^3)))))/3
    p2 = 2*p*cos(a) - L2/(3*L3)
    return(cbind(p1 = p2 + rd, p2))
  }

  wrd = function(rd, MaxIter=50) {
    pw = n1*n2/(n1 + n2)
    for (i in 1:MaxIter) {
      p = p1p2rd(rd)
      w = pw/sum(pw)
      r1 = sum(w*p[, 1])
      r2 = sum(w*p[, 2])
      w = 1/(r1*(1 - r1)/(r2*(1 - r2))/n1 + 1/n2)
      if (sum(abs(w - pw)) < 1e-8) break
      pw = w
    }
    return(w)
  }

  w = n1*n2/(n1 + n2) # MH weight for initial guess
  pRD = sum(w/sum(w)*(p1 - p2)) # initial guess with MH weight
  for (i in 1:50) {
    w2 = wrd(pRD)
    fRD = sum(w2/sum(w2)*(p1 - p2))
    if (abs(fRD - pRD) < eps) break
    pRD = fRD
  }

  Obj = function(rd) {
    p = p1p2rd(rd)
    p1 = p[, 1]
    p2 = p[, 2]
    v = (p1*(1 - p1)/n1 + p2*(1 - p2)/n2)*(n1 + n2)/(n1 + n2 - 1)
    v[v < 0] = 0
    w = wrd(rd)
    return(sum(w*(rd - fRD))^2/sum(w*w*v) - v0) # for uniroot
  }

  options(warn=-1)
  if (fRD < -1 + eps) { lower = -1
  } else { lower = uniroot(Obj, c(max(-1, fRD - 1) + eps, fRD - eps))$root }
  if (fRD > 1 - eps) { upper = 1
  } else { upper = uniroot(Obj, c(fRD + eps, 1 - eps))$root }
  options(warn=0)

  for (i in 2:nr) {
    Res1 = rbind(Res1, RDmn1(y1[i], n1[i], y2[i], n2[i], conf.level=conf.level))
  }
  rownames(Res1) = rownames(d0)    
  p1s = sum(w2/sum(w2)*p1)
  p2s = sum(w2/sum(w2)*p2)
  Res2 = c(p1=p1s, p2=p2s, RD=fRD, lower=lower, upper=upper)
  Res1 = list(Strata=Res1, Common=Res2)
  return(Res1)
}
