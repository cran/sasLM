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
  Tval[is.na(r1$DFr2) | is.nan(Tval)] = NA
  Pval[is.na(r1$DFr2) | is.nan(Pval)] = NA

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
  x = ModelMatrix(Formula, Data)
  if (!attr(terms(x), "response")) stop("Dependent variable should be provided")

  mf = model.frame(Formula, Data)
  if (!is.numeric(mf[,1])) stop("Dependent variable should be numeric!")

  cn0 = colnames(mf)
  nc0 = ncol(mf)
  vc0 = rep(FALSE, nc0)        # vector continuous 0
  for (i in 1:nc0) vc0[i] = is.numeric(mf[,i])
  cv0 = cn0[vc0][-1]           # continuous x variable names
  ncv0 = length(cv0)

  Labels = labels(terms(x))
  nLabel = length(Labels)

  tN = attr(terms(x), "factors") # table for nest information

  fIntercept = attr(x$terms, "intercept") # intercept of original formula
  if (fIntercept) { # if original formula does not have intercept, add intercept
    X = x$X
  } else {
    X = cbind(1, x$X)
    colnames(X) = c("(Intercept)", colnames(x$X))
    x$X = X
    x$assign = c(0, x$assign)
    for (i in 1:nLabel) x$termIndices[[i]] = x$termIndices[[i]] + 1
    attr(x$terms, "intercept") = 1
  }

  XpX = crossprod(X)
  nc = ncol(XpX)

  L0 = matrix(NA, nrow=nc, ncol=nc)
  dimnames(L0) = dimnames(XpX)

  L0[,1] = 1 # intercept column
  L1 = as.numeric(XpX[1,] > 0) # intercept row, columns with observation count > 0

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

  if (!fIntercept) { # if original formula does not have intercept
    L0 = L0[-1, -1]
    attr(L0, "nonEst") = iNE[-1]
  }
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

GrpCode = function(NamesInOrder, rPDIFF, conf.level=0.95)
{
  nLevel =  length(NamesInOrder)
  nL = nrow(rPDIFF)

  rowNames = rownames(rPDIFF)
  Names = strsplit(rowNames, " - ")

  m1 = matrix(NA, nrow=nLevel, ncol=nLevel)
  colnames(m1) = NamesInOrder
  rownames(m1) = NamesInOrder

  iPr = ncol(rPDIFF)
  for (k in 1:nL) {
    cX = Names[[k]][1]
    cY = Names[[k]][2]
    if (rPDIFF[k, iPr] < (1 - conf.level)) {
      m1[cX, cY] = 0
      m1[cY, cX] = 0
    }
  }

  cI = 1
  for (j in 1:nLevel) {
    fUsed = FALSE
    for (i in j:nLevel) {
      if (is.na(m1[i, j])) {
        m1[i,j] = cI
        fUsed = TRUE
      } else {
        break
      }
    }
    if (j > 1) {
      if (m1[i, j - 1] > 0) m1[j:nLevel, j] = m1[i, j - 1]
    }
    if (fUsed) cI = cI + 1
  }

  m1[is.na(m1)] = 0
  for (j in 2:nLevel) {
    fInc = TRUE
    for (i in j:nLevel) {
      if (m1[i,j] > 0 & m1[i, j - 1] == 0) fInc = FALSE
    }
    if (fInc) {
      for (i in j:nLevel) if (m1[i,j] > 0) m1[i,j] = m1[i, j - 1]
    } else {
      for (i in j:nLevel) if (m1[i,j] > 0) m1[i,j] = m1[j - 1, j - 1] + 1      
    }
  }

  gCode = rep("", nLevel)
  for (i in 1:nLevel) {
    tv = sort(unique(m1[i,]))
    t1 = paste(rep(" ", min(tv[tv > 0]) - 1), collapse="")
    t2 = paste(unique(LETTERS[m1[i,1:i]]), collapse="")
    t3 = paste(rep(" ", max(m1) - nchar(t1) - nchar(t2)), collapse="")
    gCode[i] = paste0(t1, t2, t3)
  }
  return(gCode)
}

plotDiff = function(lSM, m0, conf.level=0.95, ...)
{
  nL = nrow(m0)
  nc = ncol(m0)

  nLevel =  length(lSM)

  rowNames = rownames(m0)
  Names = strsplit(rowNames, " - ")
  Lnames = sort(unique(unlist(Names)))
  
  hDL0 = max((m0[,3] - m0[,2])/4) # Half DL is necessary
  if (is.na(hDL0) | is.nan(hDL0)) stop("SE is not available!")
  xmin = min(lSM) - hDL0
  xmax = max(lSM) + hDL0

  plot(0, 0, xlim=c(xmin, xmax), ylim=c(xmin, xmax), xlab="", ylab="", type="n", ...)
  abline(a=0, b=1, lty=2)
  abline(h=lSM, lty=3)
  abline(v=lSM, lty=3)
  text(x=xmax, y=lSM, labels=names(lSM))
  text(x=lSM, y=xmin, labels=names(lSM))

  m2Col = c("x", "y", "hDL")
  m2 = matrix(nrow=nL, ncol=length(m2Col))
  colnames(m2) = m2Col

  for (i in 1:nL) {
    m2[i,"x"] = lSM[Names[[i]][1]][[1]]
    m2[i,"y"] = lSM[Names[[i]][2]][[1]]
  }
  m2[,"hDL"] = (m0[,3] - m0[,2])/4  # PE -> PE/sqrt(2) -> DL/PE/sqrt(2) -> DL/PE/sqrt(2)/sqrt(2)*PE -> DL/2

  for (i in 1:nL) {
    if (m0[i, 1] > 0) { # For left upper region plot, PE of difference should be negative
      temp = m2[i, "x"]
      m2[i, "x"] = m2[i, "y"]
      m2[i, "y"] = temp
    }
    hDL = m2[i,"hDL"] # Half DL, -> DL/sqrt(2) if 45 degree rotation
    if (m0[i, nc] < 1 - conf.level) { Col="blue"
    } else { Col="red" }

    points(m2[i,"x"], m2[i,"y"], pch=16, col=Col)
    lines(c(m2[i,"x"] - hDL, m2[i,"x"] + hDL), c(m2[i,"y"] + hDL, m2[i,"y"] - hDL), col=Col, lwd=2)
  }
}

plotDunnett = function(m0, ...)
{
  nL = nrow(m0)
  rowNames = rownames(m0)
  Names = strsplit(rowNames, " - ")
  xNames = vector(length=nL)
  for (i in 1:nL) xNames[i] = Names[[i]][1]

  pArg = names(list(...))
  if ("xlab" %in% pArg) { Xlab = list(...)$xlab
  } else { Xlab = "Levels" }
  if ("ylab" %in% pArg) { Ylab = list(...)$ylab
  } else { Ylab = "Difference from the reference level" }
  if ("main" %in% pArg) { Title = list(...)$main
  } else { Title = paste("Difference from the reference level:", Names[[1]][2]) }
  
  ymax = max(abs(m0[,2:3]))
  
  plot(x=0, y=0, xlim=c(0.5, nL+0.5), ylim=c(-ymax, ymax), type="n", xaxt="n", xlab=Xlab, ylab=Ylab, ...)
  abline(h=0)
  axis(1, at=1:nL, labels=xNames)
  points(x=1:nL, y=m0[,1], pch=16)
  for (i in 1:nL) arrows(x0=i, y0=m0[i,2], x1=i, y1=m0[i,3], length=0.1, angle=90, code=3)
}

CheckAlias = function(Formula, Data)
{
	Res = list(Model = Formula)
  QR = qr(model.matrix(Formula, Data)) 
	Rank = QR$rank
	np = dim(QR$qr)[2]

	if (is.null(np) || Rank == np) {
	  Aliased = NULL
	} else {
		ip = 1:Rank
		dn = dimnames(QR$qr)[[2]]
		Aliased = backsolve(QR$qr[ip, ip], QR$qr[ip, -ip, drop=F])
		dimnames(Aliased) = list(dn[ip], dn[-ip])
		Aliased = t(zapsmall(Aliased))
	}
	Res$Complete = Aliased

  return(Res)
}

