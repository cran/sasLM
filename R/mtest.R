mtest = function(m1, s1, n1, m0, s0, n0, conf.level=0.95)
{
  se1 = s1/sqrt(n1)
  se0 = s0/sqrt(n0)
  ci1 = m1 + c(-1, 1)*qt(0.5 + conf.level/2, n1 - 1)*se1
  ci0 = m0 + c(-1, 1)*qt(0.5 + conf.level/2, n0 - 1)*se0
  T1 = rbind(c(m1, s1, n1, se1, ci1), c(m0, s0, n0, se0, ci0))
  colnames(T1) = c("Mean", "SD", "N", "SEM", "LCL", "UCL")
  rownames(T1) = c("Test", "Control")

  r1 = tmtest(m1, s1, n1, m0, s0, n0, conf.level=conf.level, var.equal=T)
  r2 = tmtest(m1, s1, n1, m0, s0, n0, conf.level=conf.level, var.equal=F)

  PE = rep(m1 - m0, 2)
  SE = c(r1$stderr, r2$stderr)
  LCL = c(r1$conf.int[1], r2$conf.int[1])
  UCL = c(r1$conf.int[2], r2$conf.int[2])
  Df = c(r1$parameter, r2$parameter)
  tval = c(r1$statistic, r2$statistic)
  pval = c(r1$p.value, r2$p.value)
  T2 = cbind(PE, SE, Df, LCL, UCL, tval, pval)
  colnames(T2) = c("PE", "SE", "Df", "LCL", "UCL", "t value", "Pr(>|t|)")
  rownames(T2) = c("If Equal Var.", "If Unequal Var.")
  attr(T2, "heading") = "Choose one row according to the variance equality assumption."
  class(T2) = "anova"

  Variance = c(s1^2, s0^2)
  VarLCL = c((n1 - 1)*s1^2/qchisq(0.5 + conf.level/2, n1 - 1), (n0 - 1)*s0^2/qchisq(0.5 + conf.level/2, n0 - 1))
  VarUCL = c((n1 - 1)*s1^2/qchisq(0.5 - conf.level/2, n1 - 1), (n0 - 1)*s0^2/qchisq(0.5 - conf.level/2, n0 - 1))
  T4 = cbind(Variance, VarLCL, VarUCL)
  rownames(T4) = c("Test", "Control")

  r3 = vtest(s1^2, n1, s0^2, n0, conf.level=conf.level)
  T3 = cbind(r3$estimate, r3$conf.int[1], r3$conf.int[2], r3$parameter[1], r3$parameter[2], r3$statistic, r3$p.value)
  colnames(T3) = c("PE", "LCL", "UCL", "Num Df", "Denom Df", "F value", "Pr(>F)")
  rownames(T3) = "Var. Ratio"
  class(T3) = "anova"

  Res = list(T1, T2, T4, T3)
  names(Res) = c("Statistics by group", "Two groups t-test for the difference of means", "Variance confidence intervals by group", "F-test for the ratio of variances")
  return(Res)
}

