est = function(L, X, rx, conf.level=0.95)
{
  PE = L %*% rx$coefficients
  
  if (rx$DFr > 0) {
    Var = L %*% rx$g2 %*% t(L) * rx$SSE/rx$DFr
    SE = sqrt(diag(Var))
    Tval = PE/SE
    Pval = 2*(1 - pt(abs(Tval), rx$DFr))
    DL = qt(0.5 + conf.level/2, rx$DFr)*SE
    LL =  PE - DL
    UL =  PE + DL
  } else {
    SE = NA
    Tval = NA
    Pval = NA
    LL = NA
    UL = NA
  }
  
  Res = cbind(PE, LL, UL, SE, Tval, rx$DFr, Pval)
  colnames(Res) = c("Estimate", "Lower CL", "Upper CL", "Std. Error", "t value", "Df", "Pr(>|t|)")
  attr(Res, "Estimability") = estmb(L, X, rx$g2)
  return(Res)
}

