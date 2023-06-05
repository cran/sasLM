ztest = function(m1, s1, n1, m0, s0, n0, conf.level=0.95, nullHypo=0)
{
  sampMeans = c(m1, m0)
  PE = m1 - m0
  SE = sqrt(s1^2/n1 + s0^2/n0)
  zval = (PE - nullHypo)/SE
  pval = 2*pnorm(-abs(zval))
  ci = PE + c(-1, 1)*qnorm(0.5 + conf.level/2)*SE
  attr(zval, "names") = "z"
  attr(ci, "conf.level") = conf.level
  attr(sampMeans, "names") = c("Test", "Control")
  attr(nullHypo, "names") = "difference in means"

  Res = list()
  Res$statistic = zval
  Res$p.value = pval
  Res$conf.int = ci
  Res$estimate = sampMeans
  Res$null.value = nullHypo
  Res$stderr = SE
  Res$alternative = "two.sided"
  Res$method = "Two Groups Z test"
  Res$data.name = "mean of the test (m1) and mean of the control (m0)"
  class(Res) = "htest"
  return(Res)
}
