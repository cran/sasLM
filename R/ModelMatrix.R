ModelMatrix = function(Formula, Data, NOINT=FALSE, KeepOrder=FALSE)
{
  OldOpt = options(contrasts=c("contr.treatment", "contr.poly"))
  on.exit(options(OldOpt))
  mf = model.frame(Formula, Data)
  Terms = terms(Formula, keep.order=KeepOrder)
  Labels = attr(Terms, "term.labels")

  if (!NOINT) {
    L = matrix(rep(1, nrow(mf)), ncol=1)
    colnames(L) = "(Intercept)"
    vAssign = 0
    termIndices = list()
    termIndices[[1]] = 1
  } else {
    L = NULL
    vAssign = NULL
    termIndices = list()
  }

  nLabel = length(Labels)
  for (i in 1:nLabel) {
    Lp = model.matrix(formula(paste("~", Labels[i], "- 1")), mf)
    vAssign = c(vAssign, rep(i, ncol(Lp)))
    termIndices[[i + !NOINT]] = (1:length(vAssign))[vAssign==i]
    tCol = colnames(L) # In case of continuous variable, colname does not carry.
    ColName2 = sort(colnames(Lp))
    L = cbind(L, Lp[, ColName2])
    colnames(L) = c(tCol, ColName2)
  }

  if (!NOINT) { 
    names(termIndices) = c("(Intercept)", Labels)
  } else { 
    names(termIndices) = Labels 
  }
  return(list(X=as.matrix(L), terms=Terms, termIndices=termIndices, assign=vAssign))
}
