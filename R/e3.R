e3 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, NOINT=FALSE, KeepOrder=FALSE)
  nc = NCOL(x$X)
  XpX = crossprod(x$X)
  M = getM(XpX)
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

    M1 = M
    if (!is.null(Col2)) {
      R1 = t(M[Col1, , drop=FALSE])
      R2 = t(M[Col2, , drop=FALSE])
      if (sum(abs(R1[Col1,])) > eps & sum(abs(R2[Col2,])) > eps) {
        R1 = R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1)
        M1[Col1,] = t(R1)
      }
    } else {
      M1 = pivotJ(M1, Col1, clear=FALSE)
    }
    M1 = pivotJ(M1, Col0, clear=TRUE)
    R1 = t(M1[Col1, , drop=FALSE])
    L = cbind(L, R1)
  }
  M[1:nc,] = 0 # clear all
  M[1,1] = 1 # Intercept
  if (!is.null(L)) M[colnames(L),] = t(L)
  return(M)
}
