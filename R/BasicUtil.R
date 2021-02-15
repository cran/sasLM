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

sortColName = function(ColNames)
{
  sCol = strsplit(ColNames, ":")
  nc = length(sCol[[1]])
  if (nc > 1) {
    nr = length(sCol)
    tM = matrix(unlist(sCol), nrow=nr, byrow=TRUE)
    sN = "order(nchar(tM[,1]), tM[,1]"
    for (i in 2:nc) sN = paste0(sN, ", nchar(tM[,",i,"]), tM[,", i,"]")
    sN = paste0(sN, ")")
    Ord = eval(parse(text=sN))
    Res = ColNames[Ord]
  } else {
    Res = ColNames
  }
  return(Res)
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
#  attr(Parameter, "R2") = r1$R2
#  attr(Parameter, "R2ADJ") = r1$R2ADJ
  return(Parameter)
}

lsm0 = function(x, rx, Formula, Data, conf.level=0.95)
{
  Lx = llsm0(Formula, Data)

  b = rx$coefficients[colnames(Lx)]
  PE = Lx %*% b
  if (rx$DFr > 0) {
    Var = Lx %*% rx$g2 %*% t(Lx) * rx$SSE/rx$DFr
    SE = sqrt(diag(Var))
  } else {
    SE = NA
  }

  DE = qt((1 + conf.level)/2, rx$DFr)*SE 
  LL = PE - DE
  UL = PE + DE

  Res = cbind(PE, LL, UL, SE, rx$DFr)
  colnames(Res) = c("LSmean", "LowerCL", "UpperCL", "SE", "Df")
  return(Res)
}

llsm0 = function(Formula, Data)
{
  x = ModelMatrix(Formula, Data)
  nc = NCOL(x$X)
  XpX = crossprod(x$X)

  mf = model.frame(Formula, Data)

  Labels = labels(terms(x))
  nLabel = length(Labels)
  LLabel = strsplit(Labels, ":")

  RowNames = rownames(XpX)
  ine = which(diag(XpX) == 0) # index of not estimable

  Lx = matrix(NA, nrow=nc, ncol=nc)
  dimnames(Lx) = dimnames(XpX)

  vOrd = colSums(attr(x$terms, "factors") > 0) # orders of variable
  vCont = rep(FALSE, nLabel)                   # continuous ?

## First column (intercept)
  Lx[,1] = 1

## Block Diagnoal, check continuous variable, fill first row
  for (i in 1:nLabel) {
    cI = x$termIndices[Labels[i]][[1]]
    if (length(cI) == 1) {
      for (j in i:nLabel) if (LLabel[i][[1]] %in% LLabel[j][[1]]) vCont[j] = TRUE
    } else {
      Lx[cI, cI] = diag(nrow=length(cI)) # Self is always 1.
    }
  ## First row (intercept)
    cnLevel = length(setdiff(cI, ine))
    Lx[1, cI] = 1/cnLevel
    Lx[1, ine] = 0
  }
  
## Contained or Containing
  for (i in 2:nc) {
    for (j in 2:nc) {
      if (i == j | !is.na(Lx[i,j])) next # Only for non-diagonal
      cRow = strsplit(RowNames[i], ":")[[1]]
      if (length(cRow) > 2) {
        Lx[i,] = 0 # does not support
        Lx[i,1] = NA
      }
      cCol = strsplit(RowNames[j], ":")[[1]]

      cIs = x$termIndices[Labels[x$assign[i]]][[1]]
      cJs = x$termIndices[Labels[x$assign[j]]][[1]]

      if (all(cRow %in% cCol)) {
        Lx[i, j] = length(cCol)    # Containing Flag, will be changed later
        Lx[setdiff(cIs, i), j] = 0 # This should be filled with zero
      } else if (all(cCol %in% cRow)) {
        Lx[i,j] = 1                # Contained is always 1.
        Lx[i, setdiff(cJs, j)] = 0 # This should be filled with zero
      } else if (any(cRow %in% cCol)) {
        cap = intersect(cRow, cCol)
        lcn = strsplit(colnames(Lx[,cJs]), ":")
        vcJ = rep(FALSE, length(cJs))
        for (k in 1:length(lcn)) if (all(cap %in% lcn[[k]])) vcJ[k] = TRUE
        Lx[i, cJs[vcJ]] = 1/sum(vcJ)
        Lx[i, cJs[!vcJ]] = 0
      }
    }
  }

## Containing Flag (>1) to 1/(j*k)
  for (i in 2:nc) {
    for (j in 1:nLabel) {
      cI = x$termIndices[Labels[j]][[1]]
      tI = which(Lx[i,] > 1)
      cJ = setdiff(intersect(cI, tI), ine)
      Lx[i, cJ] = 1/length(cJ)
    }
  }

## Remaining values from the intercept row
  for (j in 2:nc) {
    iL = x$assign[j]
    if (length(LLabel[iL][[1]]) == 1) {
      Lx[is.na(Lx[,j]), j] = Lx[1, j]
    } else {
      tL = strsplit(RowNames[j], ":")[[1]]
      nL = length(tL)

      oN = paste(tL[-nL], collapse=":")
      oN2 = paste(tL[-1], collapse=":")
      oN3 = paste(tL[-2], collapse=":")
      if (oN %in% RowNames) {
        for (i in which(is.na(Lx[,j]))) Lx[i, j] = Lx[oN, j]*Lx[i, oN]
      } else if (oN2 %in% RowNames) {
        for (i in which(is.na(Lx[,j]))) Lx[i, j] = Lx[oN2, j]*Lx[i, oN2]
      } else if (oN3 %in% RowNames) {
        for (i in which(is.na(Lx[,j]))) Lx[i, j] = Lx[oN3, j]*Lx[i, oN3]
      }
    }
  }

## Handle continuous (one level) variables and its interactions
  for (i in 1:nLabel) {
    cI = x$termIndices[Labels[i]][[1]]
    tL = strsplit(Labels[i], ":")[[1]]
    if (length(cI) == 1 &  length(tL) == 1) {
      Lx[,cI] = mean(mf[,Labels[i]])
    }
    if (i %in% which(vCont) & length(tL) > 1) {
      cName = colnames(Lx[,cI])
      for (j in 1:length(cName)) {
        tL2 = strsplit(cName[j], ":")[[1]]
        cVec = Lx[, tL2[1]]
        for (k in 2:length(tL2)) cVec = cVec * Lx[, tL2[k]]
        Lx[,cName[j]] = cVec
      }
    }
  }

  Lx[,ine] = 0
  Lx[ine, ] = NaN
  return(Lx[!is.na(Lx[,1]),])
}
