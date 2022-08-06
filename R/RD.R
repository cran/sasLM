RD = function(y1, n1, y2, n2, conf.level=0.95) # Risk Difference
{
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | any(c(n1, n2) == 0)) stop("Check the input!")
  p1 = y1/n1
  p2 = y2/n2
  pe = p1 - p2
  se = sqrt(p1*(1 - p1)/n1 + p2*(1 - p2)/n2)            # SE of pe
  z.crit = qnorm(0.5 + conf.level/2)
  lower = pe - z.crit*se
  upper = pe + z.crit*se
  Res = data.frame(p1=p1, p2=p2, RD=pe, SE=se, lower=lower, upper=upper)
  return(Res)
}
