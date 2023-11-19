seqCI = function(bi, ti, Zval, conf.level=0.95)
{
  bi[length(bi)] = Zval
  target = 0.5 + c(-1, 1)*conf.level/2
  LL = Drift(bi, ti, target[1])
  UL = Drift(bi, ti, target[2])
  return(c(lower=LL, upper=UL))
}
