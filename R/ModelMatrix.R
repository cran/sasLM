ModelMatrix = function(Formula, Data, KeepOrder=FALSE)
{
  if ("(Intercept)" %in% colnames(model.matrix(Formula, Data))) { 
    fIntercept = 1
  } else { 
    fIntercept = 0 
  }
  
  OldOpt = options(contrasts=c("contr.treatment", "contr.poly"))
  on.exit(options(OldOpt))
  mf = model.frame(Formula, Data)
  Terms = terms(Formula, data=Data, keep.order=KeepOrder)
  Labels = attr(Terms, "term.labels")

  if (fIntercept) {
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
    termIndices[[i + fIntercept]] = (1:length(vAssign))[vAssign==i]
    tCol = colnames(L) # In case of continuous variable, colname does not carry.
    ColName2 = sortColName(colnames(Lp))
    L = cbind(L, Lp[, ColName2])
    colnames(L) = c(tCol, ColName2)
  }

  if (fIntercept) { 
    names(termIndices) = c("(Intercept)", Labels)
  } else { 
    names(termIndices) = Labels 
  }
  return(list(X=as.matrix(L), terms=Terms, termIndices=termIndices, assign=vAssign))
}
