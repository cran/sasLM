ScoreCI = function(y, n, conf.level=0.95)
{
  if (any(c(y, n, n - y) < 0) | any(n == 0)) stop("Check the input!")
  z = qnorm(0.5 + conf.level/2)
  z2 = z^2
  mp = (y + 0.5*z2)/(n + z2) # mid point but not point estimate
  eb = z/(n + z2)*sqrt(y*(n - y)/n + z2/4) # error bound
  Res = data.frame(Prop = y/n, lower = mp - eb, upper = mp + eb)
  return(Res)
}
