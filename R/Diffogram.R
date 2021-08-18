Diffogram = function(Formula, Data, Term, conf.level=0.95, adj="lsd", ...)
{
  if (!attr(terms(Formula, data = Data), "response"))  stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")
  rx = lfit(x, y)

  if (missing(Term) & length(labels(terms(x))) == 1) Term = labels(terms(x))

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

  r1 = est(L1, x$X, rx, conf.level=conf.level, adj="lsd") # for LSMeans, do not adjust.
  r1 = r1[order(r1[,1], decreasing=TRUE),]

  m0 = PDIFF(Formula, Data, Term, conf.level=conf.level, adj=adj)
  if (tolower(adj) == "dunnett") {
    plotDunnett(m0, ...)
  } else {
    plotDiff(r1[,1], m0, conf.level=conf.level, ...)
  }
}