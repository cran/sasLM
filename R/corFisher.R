corFisher = function(x, y, conf.level=0.95, rho=0)
{ # Fisher RA. Statistical Methods for Research Workers. 14e. 1973
  n = length(x)
  r1 = cor(x, y)    # sample correlation
  z1 = atanh(r1)    # = 0.5*log((1 + r1)/(1 - r1))
  b0 = r1/(n - 1)/2 # bias to correct

  zhat = z1 - b0
  zcut = qnorm(0.5 + conf.level/2)
  seZ = 1/sqrt(n - 3)

  rhat = tanh(zhat)
  lower = tanh(zhat - zcut*seZ)
  upper = tanh(zhat + zcut*seZ)

  zn = z1 - atanh(rho) - rho/(n - 1)/2 # z under null hypothesis
  p = 2*pnorm(-abs(zn), mean=0, sd=seZ)

  Res = data.frame(N=n, r=r1, Fisher.z=z1, bias=b0, rho.hat=rhat, 
                   lower=lower, upper=upper, rho0=rho, p.value=p)
  attr(Res, "conf.level") = conf.level
  return(Res)
}

