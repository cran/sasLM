e3 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, NOINT=FALSE, KeepOrder=FALSE)
  nc = NCOL(x$X)
  XpX = crossprod(x$X)
  M = G2SWEEP(XpX, Augmented=FALSE, eps=eps) %*% XpX  
  rownames(M) = paste0("L", 1:ncol(XpX))
  Labels = labels(terms(x))
  nLabel = length(Labels)
  LLabel = strsplit(Labels, ":")

  L = NULL
  for (i in 1:nLabel) {
    Label1 = Labels[i]
    Label2 = NULL
    for (j in 1:nLabel) {
      if (i != j & all(LLabel[[i]] %in% LLabel[[j]])) Label2 = c(Label2, Labels[j])
    }

    Col1 = x$termIndices[Label1][[1]]
    Col2 = NULL
    if (length(Label2) > 0) {
      for (j in 1:length(Label2)) Col2 = c(Col2, x$termIndices[Label2[j]][[1]])
    }
    Col0 = setdiff(1:nc, c(Col1, Col2))

    if (!is.null(Col2)) {
      Ms = pivotJ(M[c(Col1, Col2), , drop=FALSE], Col0, clear=TRUE)
      Ms = pivotJ(Ms, Col1, clear=FALSE)
      bazr = which(apply(abs(Ms[, Col1, drop=FALSE]), 1, max) < eps)
      bnazr = which(apply(abs(Ms[, Col1, drop=FALSE]), 1, max) >= eps)
      bnazr = setdiff(1:length(c(Col1, Col2)), bazr)
      if (length(bnazr) == 0) {
        tL = NULL
      } else if (length(bnazr) > length(Col1)) {
        R1 = t(M[Col1, , drop=FALSE])
        R2 = t(M[Col2, , drop=FALSE])
        tL = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
      } else {
        R1 = t(Ms[bnazr, , drop=FALSE])
        R2 = t(Ms[bazr, , drop=FALSE])
        tL = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
        rownames(tL) = (rownames(M)[Col1])[1:NROW(tL)]
      }
    } else {
      tL = pivotJ(M[Col1, , drop=FALSE], Col1, clear=FALSE)
      tL = pivotJ(tL, Col0, clear=TRUE)
    }
    if (!is.null(tL)) L = rbind(L, tL)
  }
  M[1:nc,] = 0 # clear all
  M[1,1] = 1 # Intercept
  if (!is.null(L)) M[rownames(L),] = L
  return(M)
}
