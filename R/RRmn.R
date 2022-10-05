RRmn = function(d0, conf.level=0.95, eps=1e-8)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")
  v0 = qchisq(conf.level, 1)
  p1 = y1/n1
  p2 = y2/n2

  p1p2rr = function(rr) {
    A = (n1 + n2)*rr                        # eq 12
    B = n1*rr + y1 + n2 + y2*rr
    C1 = y1 + y2
    p2t = (B - sqrt(B*B - 4*A*C1))/(2*A)
    p1t = p2t*rr
    return(cbind(p1 = p2t*rr, p2 = p2t))
  }

  wrr = function(rr, MaxIter=50) {
    pw = n1*n2/(n1 + n2)
    for (i in 1:MaxIter) {
      p = p1p2rr(rr)
      w = pw/sum(pw)
      r1s = sum(w*p[, 1])
      r2s = sum(w*p[, 2])
      w = 1/((1 - r1s)/(1 - r2s)/n1 + pRR/n2)
      if (sum(abs(w - pw)) < 1e-8) break
      pw = w
    }
    return(w)
  }

  w = n1*n2/(n1 + n2) # MH weight for initial guess
  pRR = sum(w/sum(w)*p1/p2) # initial guess with MH weight
  for (i in 1:50) {
    w2 = wrr(pRR)
    RR = sum(w2/sum(w2)*p1/p2)
    if (abs(RR - pRR) < eps) break
    pRR = RR
  }

  Obj = function(rr) {
    p = p1p2rr(rr)
    r1 = p[, 1]
    r2 = p[, 2]
    w = wrr(rr)
    w = w/sum(w)
    r1s = sum(w*r1)
    r2s = sum(w*r2)
    v = (r1*(1 - r1)/n1 + RR*RR*r2*(1 - r2)/n2)*(n1 + n2)/(n1 + n2 - 1)
#    return(((r1s - r2s*RR)^2/sum(w*w*v) - v0)^2) # for nlminb
    return((r1s - r2s*RR)^2/sum(w*w*v) - v0) # for uniroot
  }

  options(warn=-1)
#  if (RR < eps) { lower = 0
#  } else { lower = nlminb(RR - eps, Obj, lower=0, upper=RR)$par }
#  upper = nlminb(RR + eps, Obj, lower=RR)$par
  if (RR < eps) { lower = 0
  } else { lower = uniroot(Obj, interval=c(eps, RR - eps))$root }
  upper = uniroot(Obj, interval=c(RR + eps, 1e9))$root
  options(warn=1)

  nr = nrow(d0)
  Res1 = RRmn1(y1[1], n1[1], y2[1], n2[1], conf.level=conf.level)
  rownames(Res1) = NULL
  if (nr > 1) {
    for (i in 2:nr) {
      Res1 = rbind(Res1, RRmn1(y1[i], n1[i], y2[i], n2[i], conf.level=conf.level))
    }
    rownames(Res1) = rownames(d0)
    p1s = sum(w2/sum(w2)*p1)
    p2s = sum(w2/sum(w2)*p2)
    Res2 = c(p1=p1s, p2=p2s, RR=RR, lower=lower, upper=upper)
    Res1 = list(Strata=Res1, Common=Res2)
  }
  return(Res1)
}
