ORcmh = function(d0, conf.level=0.95)
{
  y1 = d0[, "y1"]
  n1 = d0[, "n1"]
  y2 = d0[, "y2"]
  n2 = d0[, "n2"]
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")

  r1 = OR(y1, n1, y2, n2, conf.level=conf.level)
  
  ai = y1
  bi = n1 - y1
  ci = y2
  di = n2 - y2
  ni = n1 + n2

  pe = sum(ai*di/ni)/sum(bi*ci/ni)
  Ri = ai*di/ni
  Si = bi*ci/ni
  Pi = (ai + di)/ni
  Qi = (bi + ci)/ni
  Rsum = sum(Ri)
  Ssum = sum(Si)
  selog = sqrt(sum(Pi*Ri/(2*Rsum^2) + (Pi*Si + Qi*Ri)/(2*Rsum*Ssum) + Qi*Si/(2*Ssum^2)))

  z.crit = qnorm(0.5 + conf.level/2)
  lower = exp(log(pe) - z.crit*selog)
  upper = exp(log(pe) + z.crit*selog)
  r2 = data.frame(OR = pe, SElog = selog, lower = lower, upper = upper)
  return(list(ORs=r1, Common=r2))
}
