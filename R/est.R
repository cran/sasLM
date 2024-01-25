est = function(L, X, rx, conf.level=0.95, adj="lsd", paired=FALSE)
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
      if (paired) {
        nLevel = round(0.5 + sqrt(1 + 8*nL)/2)
        if (abs(nLevel*(nLevel - 1)/2 - nL) > 1e-8) stop("Row count of L does not match for Tukey adjustment!")
      } else {
        nLevel = nL
      }
      SE2 = SE/sqrt(2)
      Pval = 1 - ptukey(abs(PE/SE2), nLevel, rx$DFr)
      DL = qtukey(conf.level, nLevel, rx$DFr)*SE2
#    } else if (tolower(adj) == "duncan") {
#      if (paired) {
#        nLevel = round(0.5 + sqrt(1 + 8*nL)/2)
#        if (abs(nLevel*(nLevel - 1)/2 - nL) > 1e-8) stop("Row count of L does not match for Duncan adjustment!")
#      } else {
#        nLevel = nL
#      }
#      SE2 = SE/sqrt(2)
#      Pval = 1 - ptukey(abs(PE/SE2), nLevel, rx$DFr)^(1/(nLevel - 1))
#      DL = qtukey(conf.level^(nLevel - 1), nLevel, rx$DFr)*SE2
    } else if (substr(tolower(adj), 1, 3) == "bon") {
      Pval = nL*2*(1 - pt(abs(Tval), rx$DFr))
      Pval[Pval > 1] = 1
      DL = qt(1 - (1 - conf.level)/2/nL, rx$DFr)*SE
    } else if (tolower(adj) == "scheffe") {
      DFg = round(0.5 + sqrt(1 + 8*nL)/2) - 1 # Df group
      Fs = rep(0, nL)
      for (j in 1:nL) Fs[j] = cSS(L[j, , drop=F], rx)[["F value"]]
      Fval = Fs/DFg
      Pval = 1 - pf(Fval, DFg, rx$DFr)
      DL = sqrt(DFg*qf(conf.level, DFg, rx$DFr))*SE
    } else if (tolower(adj) == "dunnett") {
      ColNames = colnames(L)
      Lf = colSums(L)
      Groups = ColNames[Lf == 1]
      Control = ColNames[Lf == -nL]
      nX = diag(crossprod(X))
      nControl = nX[Control]
      nGroups = nX[Groups]
      vSD = sqrt(nGroups/(nGroups + nControl))
      mC = outer(vSD, vSD, "*")
      diag(mC) = 1
#      mC = matrix(0.5, nrow=nL, ncol=nL) + diag(0.5, nrow=nL) # for balanced data only!
#      if (exists(".Random.seed")) Saved.seed = .Random.seed
#      set.seed(5)
#        D.crit = qmvt(0.5 + conf.level/2, df=rx$DFr, corr=mC)$quantile

#      if (dim(mC)[1] < 11) {
#        D.crit = qmvt(0.5 + conf.level/2, df=rx$DFr, corr=mC, algorithm=Miwa)$quantile
#      } else {
      D.crit = qmvt(0.5 + conf.level/2, df=rx$DFr, corr=mC, seed=5)$quantile        
#      }
      Pval = rep(NA, nL)
      DL = rep(NA, nL)
      for (k in 1:nL) {
#        set.seed(5) # DescTools::DunnettTest forgot to set seed before pmvt
#        Pval[k] = 1 - pmvt(lower=rep(-abs(Tval[k]), nL), upper=rep(abs(Tval[k]), nL), df=rx$DFr, corr=mC)
#        if (dim(mC)[1] < 11) {
#          Pval[k] = 1 - pmvt(lower=rep(-abs(Tval[k]), nL), upper=rep(abs(Tval[k]), nL), df=rx$DFr, corr=mC, algorithm=Miwa)
#        } else {
        Pval[k] = 1 - pmvt(lower=rep(-abs(Tval[k]), nL), upper=rep(abs(Tval[k]), nL), df=rx$DFr, corr=mC, seed=5)
#        }
        DL[k] = D.crit*SE[k]
      }
#      if (exists("Saved.seed")) .Random.seed <<- Saved.seed
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
    Res = cbind(PE, LL, UL, SE, Fval, rx$DFr, Pval)
    colnames(Res) = c("Estimate", "Lower CL", "Upper CL", "Std. Error", "F value", "Df2", "Pr(>F)")
    attr(Res, "Estimability") = estmb(L, X, rx$g2)
  }

  return(Res)
}

