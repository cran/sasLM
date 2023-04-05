LInorm = function(x, k, conf.level=0.95, eps=1e-8)
{
  n = length(x)
  m0 = mean(x)
  s0 = sqrt(var(x)*(n - 1)/n)
  maxLL = sum(dnorm(x, mean=m0, sd=s0, log=T))

  logk = ifelse(missing(k), qf(conf.level, 1, max(n - 1, 1))/2, log(k)) # for both estimated and profile likelihood !!! DF does not depend on the number of parameters !!!
  logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05

  Obj1 = function(th) (maxLL - sum(dnorm(x, mean=th, sd=s0, log=T)) - logk)^2 # for estimated likelihood
  Obj2 = function(th) (maxLL - sum(dnorm(x, mean=m0, sd=th, log=T)) - logk)^2 # for estimated likelihood

  LLmean = nlminb(m0 - abs(m0)/10, Obj1, lower=-1e6, upper=m0)$par
  ULmean = nlminb(1.1*m0, Obj1, lower=m0, upper=1e6)$par
  
  LLsd = nlminb(s0*0.9, Obj2, lower=0, upper=s0)$par
  ULsd = nlminb(s0*1.1, Obj2, lower=s0, upper=1e6)$par
  
  Res = cbind(PE=c(m0, s0), LL=c(LLmean, LLsd), UL=c(ULmean, ULsd))
  rownames(Res) = c("mean", "sd")
  attr(Res, "k") = exp(logk)
  return(Res)
}
