tmtest = function(m1, s1, n1, m0, s0, n0, conf.level=0.95, nullHypo=0, var.equal=F)
{
  sampMeans = c(m1, m0)
  PE = m1 - m0
  if (var.equal) {
    SPsq = ((n1 - 1)*s1^2 + (n0 - 1)*s0^2)/(n1 + n0 - 2)
    SE = sqrt(SPsq*(1/n1 + 1/n0))
    Df = n1 + n0 - 2
  } else {
    ns = c(n1, n0)
    vars = c(s1^2, s0^2)/ns
    SE = sqrt(sum(vars))
    Df = sum(vars)^2/sum(vars^2/(ns - 1))
  }
  tval = (PE - nullHypo)/SE
  pval = 2*pt(-abs(tval), Df)
  ci = PE + c(-1, 1)*qt(0.5 + conf.level/2, Df)*SE
  attr(tval, "names") = "t"
  attr(Df, "names") = "df"
  attr(ci, "conf.level") = conf.level
  attr(sampMeans, "names") = c("Test", "Control")
  attr(nullHypo, "names") = "difference in means"

  Res = list()
  Res$statistic = tval
  Res$parameter = Df
  Res$p.value = pval
  Res$conf.int = ci
  Res$estimate = sampMeans
  Res$null.value = nullHypo
  Res$stderr = SE
  Res$alternative = "two.sided"
  Res$method = ifelse(var.equal, "Two Sample t-test", "Welch Two Sample t-test")
  Res$data.name = "mean of the test (m1) and mean of the control (m0)"
  class(Res) = "htest"
  return(Res)
}
