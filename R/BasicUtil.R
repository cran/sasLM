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
  if ("(Intercept)" %in% colnames(r1$g2)) {
    DF = c(r1$rank - 1, r1$DFr, nObs - 1)
  } else {
    DF = c(r1$rank, r1$DFr, nObs)
  }
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
  
  if ("(Intercept)" %in% colnames(r1$g2)) {
    rownames(ANOVA) = c("MODEL", "RESIDUALS", "CORRECTED TOTAL")
  } else {
    rownames(ANOVA) = c("MODEL", "RESIDUALS", "UNCORRECTED TOTAL")    
  }
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
  np = ncol(X)
  Est = r1$coefficients

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

  Est[is.na(r1$DFr2)] = NA
  bSE[is.na(r1$DFr2)] = NA
  Tval[is.na(r1$DFr2)] = NA
  Pval[is.na(r1$DFr2)] = NA

  ESTM = estmb(diag(np), X, r1$g2)
  if (sum(ESTM) == np) {
    Parameter = cbind(Est, bSE, r1$DFr2, Tval, Pval)
    colnames(Parameter) = c("Estimate", "Std. Error", "Df", "t value", "Pr(>|t|)")    
  } else {
    Parameter = cbind(Est, ESTM, bSE, r1$DFr2, Tval, Pval)
    colnames(Parameter) = c("Estimate", "Estimable", "Std. Error", "Df", "t value", "Pr(>|t|)")
  }
  rownames(Parameter) = colnames(X)
  class(Parameter) = "anova"
#  attr(Parameter, "R2") = r1$R2
#  attr(Parameter, "R2ADJ") = r1$R2ADJ
  return(Parameter)
}

lsm0 = function(x, rx, Formula, Data, conf.level=0.95, hideNonEst=TRUE)
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
  if (hideNonEst) Res[attr(Lx, "nonEst"),] = NA
  return(Res)
}

llsm0 = function(Formula, Data)
{
# Find continuous variables at model frame
  mf = model.frame(Formula, Data)
  cn0 = colnames(mf)
  nc0 = ncol(mf)
  vc0 = rep(FALSE, nc0)        # vector continuous 0
  for (i in 1:nc0) vc0[i] = is.numeric(mf[,i])
  cv0 = cn0[vc0][-1]           # continuous x variable names
  ncv0 = length(cv0)

  x = ModelMatrix(Formula, Data)
  Labels = labels(terms(x))
  nLabel = length(Labels)

  XpX = crossprod(x$X)
  nc = ncol(XpX)

  L0 = matrix(NA, nrow=nc, ncol=nc)
  dimnames(L0) = dimnames(XpX)

  L0[,1] = 1

  tN = attr(terms(x), "factors") # table for nest information
  L1 = as.numeric(XpX[1,] > 0)

# First row (intercept)
  for (i in 1:nLabel) {
    cI = x$termIndices[Labels[i]][[1]]
    if (sum(tN[,Labels[i]] > 1) > 0) {  # nested
      nestF = rownames(tN)[tN[,Labels[i]] > 1]
      nestF = paste(nestF, collapse=":")
      nfCI = x$termIndices[[nestF]]
      nnfCI = length(nfCI)
      nCI = length(cI)
      sG = nCI/nnfCI
      if (sG > 1) {
        for (j in 1:nnfCI) {
          ccI = cI[((j - 1)*sG + 1):(j*sG)]
          L0[1, ccI] = 1/sum(L1[ccI])*L0[1, nfCI[j]]
        }
      } else {
        L0[1, cI] = 1/sum(L1[cI])
      }
    } else {
      L0[1, cI] = 1/sum(L1[cI])
    }
  }

  for (i in 1:nLabel) {
    for (j in 1:nLabel) {
      Ls = bL(i, j, x, L0[1,,drop=FALSE])
      L0[rownames(Ls), colnames(Ls)] = Ls
    }
  }

# Fill single continuous columns with mean
  if (ncv0 > 0) {
    for (i in 1:ncv0) L0[, cv0[i]] = mean(mf[, cv0[i]])
  }

# Two or more interaction columns with continuous variable
  cn1 = colnames(XpX)
  scn1 = strsplit(cn1, ":")
  names(scn1) = cn1

  for (i in 1:nc) {
    ccn = scn1[[i]]           # splitted current column name
    nccn = length(ccn)        # number of ccn
    if (nccn == 1) next       # already treated above, no more need to consider for single name
    cap = intersect(cv0, ccn)
    ncap = length(cap)
    if (ncap == 0) next       # no intersection -> no contiuous varaible
    tv = L0[, ccn[1]]         # temporary vector
    for (j in 2:nccn) tv = tv*L0[, ccn[j]] # it must be more than 1 column (nccn >= 2)
    L0[, i] = tv
  }

# Check estimability  
  iNE = (diag(XpX) == 0) # index of not estimable
  nNE = sum(iNE)
  if (nNE > 0) {
    sNE = strsplit(cn1[iNE], ":")
    for (i in 1:nNE) {
      cStr = sNE[[i]]
      for (j in 2:nc) {
        cCol = scn1[[j]]
#        if (all(cCol %in% cStr) & length(cCol) < length(cStr) & tN[x$assign[j] + 1, x$assign[iNE][i]] == 1) iNE[j] = TRUE
        C1 = all(cCol %in% cStr) & length(cCol) < length(cStr)
        C2 = length(cCol)
        if (C2 > 1) {
          if (C1) iNE[j] = TRUE
        } else if (C2 == 1) {
           if (C1 & tN[x$assign[j] + 1, x$assign[iNE][i]] < 2) iNE[j] = TRUE
        }
      }
    }
  }

  L0[, diag(XpX) == 0] = 0 # No observation columns
  attr(L0, "nonEst") = iNE

  return(L0)
}

bL = function(m, n, x, L1)
{
  Labels = labels(terms(x))
  tLabel = Labels[m]
  cLabel = Labels[n]

  sLabel = strsplit(Labels, ":")
  fCap = intersect(sLabel[[m]], sLabel[[n]])

  XpX = crossprod(x$X)
  cn = colnames(XpX)
  iNO = (1:ncol(XpX))[diag(XpX) == 0]
  cI = x$termIndices[[tLabel]]
  cJ = setdiff(x$termIndices[[cLabel]], iNO)

  nr = length(cI)
  nc = length(cJ)
  Lx = matrix(0, nrow=nr, ncol=nc)
  dimnames(Lx) = list(cn[cI], cn[cJ])

  if (length(fCap) == 0) {
    for (i in 1:nr) Lx[i, ] = L1[1, setdiff(cJ, iNO)]
  } else if (m == n) {
    Lx = diag(1, nrow=length(cI))
    dimnames(Lx) = list(cn[cI], cn[cI])
  } else {
    tLevels = strsplit(cn[cI], ":")
    cLevels = strsplit(cn[cJ], ":")
    iLevels = intersect(unlist(tLevels), unlist(cLevels))

    tNames = rep("", nr)
    cNames = rep("", nc)
    for (i in 1:nr) tNames[i] = paste(intersect(tLevels[[i]], iLevels), collapse=":")
    for (j in 1:nc) cNames[j] = paste(intersect(cLevels[[j]], iLevels), collapse=":")
    for (i in 1:nr) for (j in 1:nc) if (tNames[i] == cNames[j]) Lx[i, j] = 1

    Lx = Lx/rowSums(Lx)
  }

  return(Lx)
}

