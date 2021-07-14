PDIFF = function(Formula, Data, Term, conf.level=0.95, adj="lsd", ref, PLOT=FALSE, reverse=FALSE)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  if (missing(Term) & length(labels(terms(x))) == 1) Term = labels(terms(x))

  rx = lfit(x, y)

  L0 = llsm0(Formula, Data)
  nc = NCOL(L0)

  Labels = labels(terms(x))
  ti = x$termIndices[Term][[1]]
  nti = length(ti)

  sTerm = strsplit(Term, ":")[[1]]
  if (length(sTerm) > 1) {
    ColNames = vector()
    for (i in 1:nti) {
      cCol = colnames(L0)[ti[i]]
      cStr = strsplit(cCol, ":")[[1]]
      tStr = vector()
      for (j in 1:length(sTerm)) {
        tStr[j] = substring(cStr[j], nchar(sTerm[j]) + 1)
      }
      ColNames[i] = paste(tStr, collapse=":")
    }
  } else {
    ColNames = substring(colnames(L0)[ti], nchar(Term) + 1)
  }

  if (nti == 0) stop(paste(Term, "term not found!"))
  if (nti == 1) return(est(L0, x$X, rx, conf.level=conf.level, adj=adj))

  nr = nti*(nti - 1)/2
  Lx = matrix(ncol=nc, nrow=nr)
  RowNames = vector(length = nr)
  iL1 = 1
  for (i in 1:nti) {
    for (j in (i + 1):nti) {
      if (j == i) next
      Lx[iL1, ] = L0[ti[i], ] - L0[ti[j], ]
      if (reverse) {
        RowNames[iL1] = paste(ColNames[j], "-", ColNames[i])        
      } else {
        RowNames[iL1] = paste(ColNames[i], "-", ColNames[j])
      }
      iL1 = iL1 + 1
    }
    if (iL1 > nr) break
  }
  if (reverse) Lx = -Lx

  if (tolower(adj) == "dunnett") { # adjust Lx
    nL = nrow(Lx)
    Names = strsplit(RowNames, " - ")
    RemoveI = NULL
    for (i in 1:nL) {
      if (ref %in% Names[[i]]) {
        if (ref == Names[[i]][1]) {
          Lx[i, ] = -1 * Lx[i,]
          RowNames[i] = paste(Names[[i]][2], "-", Names[[i]][1])
        }
      } else {
        RemoveI = c(RemoveI, i)
      }
    }
    Lx = Lx[-RemoveI,]
    RowNames = RowNames[-RemoveI]
  }

  rownames(Lx) = RowNames
  colnames(Lx) = names(rx$coefficients)

  if (tolower(adj) == "duncan") {
    PE = Lx %*% rx$coefficients
    Ms = as.vector(L0 %*% rx$coefficients)[ti]
    n = (abs(outer(rank(Ms), rank(Ms), "-")) + 1)[lower.tri(diag(length(Ms)))]

    if (rx$DFr > 0) {
      Var = Lx %*% rx$g2 %*% t(Lx) * rx$SSE/rx$DFr
      SE2 = sqrt(diag(Var))/sqrt(2)
      Tval = abs(PE)/SE2
      Pval = 1 - ptukey(Tval, n, rx$DFr)^(1/(n - 1))
      DL = qtukey(conf.level^(n - 1), n, rx$DFr)*SE2
      LL =  PE - DL
      UL =  PE + DL
    } else {
      Pval = NA
      LL = NA
      UL = NA
    }
    Res = cbind(PE, LL, UL, Pval)
    colnames(Res) = c("Estimate", "Lower CL", "Upper CL", "Pr(>|t|)")
    attr(Res, "Estimability") = estmb(Lx, x$X, rx$g2)
  } else {
    Res = est(Lx, x$X, rx, conf.level=conf.level, adj=adj)
  }
  
  class(Res) = "anova"
  
  if (PLOT) {
    L1 = L0[ti, , drop=F]
    nL = NROW(L1)
    rowNames = rownames(L1)
    newRowNames = vector(length=nL)

    sTerm = strsplit(Term, ":")[[1]]
    if (length(sTerm) > 1) {
      for (i in 1:nL) {
        cLevel0 = strsplit(rowNames[i], ":")[[1]]
        cLevel1 = vector(length=length(sTerm))
        for (j in 1:length(cLevel0)) {
          cLevel1[j] = substr(cLevel0[j], nchar(sTerm[j]) + 1, nchar(cLevel0[j]))
        }
        newRowNames[i] = paste(cLevel1, collapse=":")
      }
    } else {
      si = nchar(Term) + 1
      for (i in 1:nL) newRowNames[i] = substr(rowNames[i], si, nchar(rowNames[i]))
    }
    rownames(L1) = newRowNames

    r1 = est(L1, x$X, rx, conf.level=conf.level, adj="lsd") # for LSMeans, do not adjust.
    r1 = r1[order(r1[,1], decreasing=TRUE),]
    if (tolower(adj) == "dunnett") {
      plotDunnett(r1)
    } else {    
      Title = paste("Diffogram of", Term)
      plotDiff(r1[,1], Res, conf.level=conf.level, Title=Title)
    }
  }

  return(Res)
}
