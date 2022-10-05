RR = function(y1, n1, y2, n2, conf.level=0.95) # Relative Risk
{
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(n1*n2 == 0)) stop("Check the input!")
  p1 = y1/n1
  p2 = y2/n2
  pe = p1/p2
  selog = sqrt(1/y1 - 1/n1 + 1/y2 - 1/n2)               # SE of log(pe)
  z.crit = qnorm(0.5 + conf.level/2)
  lower = exp(log(pe) - z.crit*selog)
  upper = exp(log(pe) + z.crit*selog)
  Res = data.frame(p1=p1, p2=p2, RR=pe, SElog=selog, lower=lower, upper=upper)
  return(Res)
}
