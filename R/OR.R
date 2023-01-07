OR = function(y1, n1, y2, n2, conf.level=0.95) # Odds Ratio
{
  ta = c(y1, n1 - y1, y2, n2 - y2)
  if (any(ta < 0) | any(n1*n2 == 0)) stop("Check the input!")
  if (any(ta == 0)) {
    warning("Values are added: 0.5s to the numberator, 1s to the denominator.\n  Consider other ways too!")
    y1 = y1 + 0.5
    n1 = n1 + 1
    y2 = y2 + 0.5
    n2 = n2 + 1
  }
  o1 = y1/(n1 - y1) # odd 1
  o2 = y2/(n2 - y2) # odd 2
  pe = o1/o2
  selog = sqrt(1/y1 + 1/(n1 - y1) + 1/y2 + 1/(n2 - y2)) # SE of log(pe)
  z.crit = qnorm(0.5 + conf.level/2)
  lower = exp(log(pe) - z.crit*selog)
  upper = exp(log(pe) + z.crit*selog)
  Res = data.frame(odd1=o1, odd2=o2, OR=pe, SElog=selog, lower=lower, upper=upper)
  return(Res)  
}
