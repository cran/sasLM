est = function(L, rx)
{
  PE = L %*% rx$coefficients
  
  if (rx$DFr > 0) {
    Var = L %*% rx$g2 %*% t(L) * rx$SSE/rx$DFr
    SE = sqrt(diag(Var))
    Tval = PE/SE
    Pval = 2*(1 - pt(abs(Tval), rx$DFr))
  } else {
    SE = NA
    Tval = NA
    Pval = NA
  }
  
  Res = cbind(PE, SE, Tval, Pval)
  colnames(Res) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  class(Res) = "anova"
  return(Res)
}
