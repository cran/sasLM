LSM = function(Formula, Data, Term, conf.level=0.95, adj="lsd", hideNonEst=TRUE, PLOT=FALSE, ...)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  rx = lfit(x, y)

  if (missing(Term) & length(labels(terms(x))) == 1) Term = labels(terms(x))

  if (missing(Term)) {
    Res = lsm0(x, rx, Formula, Data, conf.level=conf.level, hideNonEst=hideNonEst)
  } else {
    L0 = llsm0(Formula, Data)
    nc = NCOL(L0)

    Labels = labels(terms(x))
    ti = x$termIndices[Term][[1]]
    nti = length(ti)

    L1 = L0[ti,,drop=F]
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

    r1 = est(L1, x$X, rx, conf.level=conf.level, adj="lsd")
    r1 = r1[order(r1[,1], decreasing=TRUE),]
    r2 = PDIFF(Formula, Data, Term=Term, conf.level=conf.level, adj=adj)

    Res = data.frame(GrpCode(rownames(r1), r2), r1[,c(1,2,3,4,6)])
    colnames(Res) = c("Group", "LSmean", "LowerCL", "UpperCL", "SE", "Df")

    if (PLOT) {
      pArg = names(list(...))
      if ("xlab" %in% pArg) { Xlab = list(...)$xlab
      } else { Xlab = Term }
      if ("ylab" %in% pArg) { Ylab = list(...)$ylab
      } else { Ylab = as.character(Formula)[2] }

      ymin = min(r1[,2])
      ymax = max(r1[,3])
      Range = ymax - ymin
      plot(0, 0, xlim=c(0.5, (nL + 0.5)), ylim=c(ymin - 0.2*Range, ymax + 0.2*Range), xlab=Xlab, ylab=Ylab, xaxt="n", type="n", ...)
      axis(1, at=1:nL, labels=rownames(r1))
      points(x=1:nL, y=r1[,1], pch=16)
      for (i in 1:nL) arrows(x0=i, y0=r1[i,2], x1=i, y1=r1[i,3], length=0.1, angle=90, code=3)
    }
  }
  return(Res)
}

