G2SWEEP = function(A, Augmented=FALSE, eps=1e-8)
{
  p = nrow(A)
  p0 = ifelse(Augmented, p - 1, p)
  r = 0
  for (k in 1:p0) {
    d = A[k,k]
#    CSS = crossprod(A[k,] - mean(A[k,]))
#    DminK = ifelse(CSS > 0, eps*CSS, eps)
    if (abs(d) < eps) { A[k,] = 0 ; A[,k] = 0 ; next }
    A[k,] = A[k,]/d
    r = r + 1
    for (i in 1:p) {
      if (i != k) {
        c0 = A[i,k] ; 
        A[i,] = A[i,] - c0*A[k,] ; 
        A[i,k] = -c0/d
      }
    }
    A[k,k] = 1/d
  }
  attr(A, "rank") = r
  return(A)
}

pivotJ = function(M, j, clear=TRUE, eps=1e-8)
{
  for (k in j) {
    if (any(abs(M[,k]) > eps) > 0) {
      Js = which(abs(as.vector(M[,k])) > eps)
      pivotRow = M[Js[1], ]/M[Js[1], k]
      nJ = length(Js)
      if (nJ > 1) for (i in 2:nJ) M[Js[i],] = M[Js[i],] - M[Js[i], k]*pivotRow
      if (clear) { 
        M[Js[1],] = 0
      } else {
        M[Js[1],] = pivotRow
      }
    }
  }
  return(M)
}

sumANOVA = function(r1, T1, SST, nObs, yName=NULL)
{
  DF = c(r1$rank - 1, r1$DFr, nObs - 1)
  SS = c(SST - r1$SSE, r1$SSE, SST)
  if (DF[2] > 0) {
    MS = c(SS[1:2]/DF[1:2], NA)
  } else {
    MS = c(SS[1]/DF[1], NA, NA)
  }
  if (MS[2] > 0 & DF[2] > 0) {
    Fval = c(MS[1]/MS[2], NA, NA)
    Pval = c(1 - pf(Fval[1], DF[1], DF[2]), NA, NA)
  } else {
    Fval = rep(NA, 3)
    Pval = rep(NA, 3)
  }

  ANOVA = cbind(DF, SS, MS, Fval, Pval)
  colnames(ANOVA) = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  rownames(ANOVA) = c("MODEL", "RESIDUALS", "CORRECTED TOTAL")
  if (!is.null(T1)) {
    rownames(T1) = paste0(" ", rownames(T1))
    ANOVA = rbind(ANOVA[1,,drop=FALSE], T1, ANOVA[2:3,])
  }
  if (!is.null(yName)) {
    attr(ANOVA, "heading") = paste("Response :", yName)
  }
  class(ANOVA) = "anova"
  return(ANOVA)
}

sumREG = function(r1, X)
{
  if (r1$DFr > 0) {
    bVar = r1$g2 %*% crossprod(X) %*% t(r1$g2) * r1$SSE/r1$DFr
    bSE = sqrt(diag(bVar))
    Tval = r1$coefficients/bSE
    Pval = 2*(1 - pt(abs(Tval), r1$DFr))
  } else {
    bSE = NA
    Tval = NA
    Pval = NA
  }
  Parameter = cbind(r1$coefficients, bSE, Tval, Pval)
  colnames(Parameter) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(Parameter) = colnames(X)
  class(Parameter) = "anova"

  return(Parameter)
}
