PDIFF = function(Formula, Data, Term, conf.level=0.95)
{
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[,1]
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
  if (nti == 1) return(est(L0, x$X, rx, conf.level=conf.level))

  nr = nti*(nti - 1)/2
  Lx = matrix(ncol=nc, nrow=nr)
  RowNames = vector(length = nr)
  iL1 = 1
  for (i in 1:nti) {
    for (j in (i + 1):nti) {
      if (j == i) next
      Lx[iL1, ] = L0[ti[i], ] - L0[ti[j], ]
      RowNames[iL1] = paste(ColNames[i],"-",ColNames[j])
      iL1 = iL1 + 1
    }
    if (iL1 > nr) break
  }

  Res = est(Lx, x$X, rx, conf.level=conf.level)
  rownames(Res) = RowNames
  printCoefmat(Res)
  invisible(Res)
}