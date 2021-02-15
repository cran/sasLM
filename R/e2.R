e2 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  X = x$X
  nc = NCOL(X)
  XpX = crossprod(X)
  M0 = G2SWEEP(XpX, Augmented=FALSE, eps=eps) %*% XpX
  rownames(M0) = paste0("L", 1:ncol(XpX))
  Labels = labels(terms(x))
  nLabel = length(Labels)
  LLabel = strsplit(Labels, ":")

  Ls = c(1, rep(0, nc - 1)) # interecept
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

    X0 = X[,Col0]
    X1 = X[,Col1]
    X2 = X[,Col2]
    Mx = X0 %*% G2SWEEP(crossprod(X0)) %*% t(X0)
    M = diag(NCOL(Mx)) - Mx

    X1pM  = crossprod(X1, M)
    X1pMX1 = X1pM %*% X1
    gX1pMX1 = G2SWEEP(X1pMX1)

    L = M0[x$termIndices[Label1][[1]],,drop=FALSE]
    L[,Col0] = 0
    L[,Col1] = gX1pMX1 %*% X1pMX1
    L[,Col2] = gX1pMX1 %*% X1pM %*% X2
    Ls = rbind(Ls, L)
  }
  rownames(Ls) = paste0("L", 1:nc)
  return(Ls)
}
