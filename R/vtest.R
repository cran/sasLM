vtest = function(v1, n1, v0, n0, ratio=1, conf.level=0.95)
{
  PE = v1/v0
  LCL = PE/qf(0.5 + conf.level/2, n1 - 1, n0 - 1)
  UCL = PE/qf(0.5 - conf.level/2, n1 - 1, n0 - 1)
  ci = c(LCL, UCL)

  Fval = v1/v0/ratio
  attr(Fval, "names") = "F"
  if (PE > 1) {
    pval = 2*(1 - pf(Fval, n1 - 1, n0 - 1))
  } else {
    pval = 2*pf(Fval, n1 - 1, n0 - 1)
  }
  Df = c(n1 - 1, n0 - 1)
  attr(PE, "names") = "ratio of variances"
  attr(Df, "names") = c("num df", "denom df")
  attr(ci, "conf.level") = conf.level
  attr(ratio, "names") = "ratio of variancess"

  Res = list()
  Res$statistic = Fval
  Res$parameter = Df
  Res$p.value = pval
  Res$conf.int = ci
  Res$estimate = PE
  Res$null.value = ratio
  Res$alternative = "two.sided"
  Res$method = "F test to compare two variances"
  Res$data.name = "variance of the test (v1) and variance of the control (v0)"
  class(Res) = "htest"
  return(Res)
}
