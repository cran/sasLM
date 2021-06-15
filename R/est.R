est = function(L, X, rx, conf.level=0.95, adj="lsd")
{
  nL = nrow(L)
  PE = L %*% rx$coefficients

  if (rx$DFr > 0) {
    Var = L %*% rx$g2 %*% t(L) * rx$SSE/rx$DFr
    SE = sqrt(diag(Var))
    Tval = PE/SE
    if (tolower(adj) == "lsd") {
      Pval = 2*(1 - pt(abs(Tval), rx$DFr))
      DL = qt(0.5 + conf.level/2, rx$DFr)*SE
    } else if (tolower(adj) == "tukey") {
      nLevel = round(0.5 + sqrt(1 + 8*nL)/2)
      if (abs(nLevel*(nLevel - 1)/2 - nL) > 1e-8) stop("Row count of L does not match for Tukey adjustment!")
      Pval = 1 - ptukey(abs(Tval*sqrt(2)), nLevel, rx$DFr)
      DL = qtukey(conf.level, nLevel, rx$DFr)*SE/sqrt(2)
    } else if (substr(tolower(adj), 1, 3) == "bon") {
      Pval = nL*2*(1 - pt(abs(Tval), rx$DFr))
      Pval[Pval > 1] = 1
      DL = qt(1 - (1 - conf.level)/2/nL, rx$DFr)*SE
    } else if (tolower(adj) == "scheffe") {
      DFg = round(0.5 + sqrt(1 + 8*NROW(L))/2) - 1 # Df group
      Fs = rep(0, nL)
      for (j in 1:nL) Fs[j] = cSS(L[j, , drop=F], rx)[["F value"]]
      Pval = 1 - pf(Fs/DFg, DFg, rx$DFr)
      DL = sqrt(DFg*qf(conf.level, DFg, rx$DFr))*SE
    } else if (tolower(adj) == "dunnett") {
      mC = matrix(0.5, nrow=nL, ncol=nL) + diag(0.5, nrow=nL)
      Pval = rep(NA, nL)
      DL = rep(NA, nL)
      for (k in 1:nL) {
        Pval[k] = 1 - pmvt(lower=rep(-abs(Tval[k]), nL), upper=rep(abs(Tval[k]), nL), df=rx$DFr, corr=mC)
        DL[k] = qmvt(0.5 + conf.level/2, df=rx$DFr, corr=mC)$quantile*SE[k]
      }
    } else {
      stop(paste0("Adjustment method ", adj, " is not supported!"))
    }
    LL =  PE - DL
    UL =  PE + DL
  } else {
    SE = NA
    Tval = NA
    Pval = NA
    LL = NA
    UL = NA
  }

  if (tolower(adj) %in% c("lsd", "tukey", "bon", "dunnett")) {
    Res = cbind(PE, LL, UL, SE, Tval, rx$DFr, Pval)
    colnames(Res) = c("Estimate", "Lower CL", "Upper CL", "Std. Error", "t value", "Df", "Pr(>|t|)")
    attr(Res, "Estimability") = estmb(L, X, rx$g2)
  } else {
    Res = cbind(PE, LL, UL, Pval)
    colnames(Res) = c("Estimate", "Lower CL", "Upper CL", "Pr(>F)")
    attr(Res, "Estimability") = estmb(L, X, rx$g2)
  }

  return(Res)
}

