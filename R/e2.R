e2 = function(x, eps=1e-8)
{
  eps2 = eps/max(abs(x$X))
  nc = NCOL(x$X)
  XpX = crossprod(x$X)
  M0 = G2SWEEP(XpX, Augmented=FALSE, eps=eps) %*% XpX
  M0[abs(M0) < eps2] = 0
  rownames(M0) = paste0("L", 1:nc)
  Labels = labels(terms(x))
  nLabel = length(Labels)
  LLabel = strsplit(Labels, ":")

  if (attr(x$terms, "intercept")) {
    re2 = c(1, rep(0, nc - 1))
  } else {
    re2 = NULL
  }

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

    X0 = x$X[, Col0]
    X1 = x$X[, Col1]
    X2 = x$X[, Col2]

    if (NCOL(X0) > 0) {
      Mx = X0 %*% G2SWEEP(crossprod(X0), Augmented=FALSE, eps=eps) %*% t(X0)
      Mx[abs(Mx) < eps] = 0
      M = diag(NCOL(Mx)) - Mx
      X1pM  = crossprod(X1, M)
    } else {
      X1pM = t(X1)
    }
    X1pMX1 = X1pM %*% X1
    gX1pMX1 = G2SWEEP(X1pMX1, Augmented=FALSE, eps=eps)

    L2 = M0[x$termIndices[Label1][[1]], , drop=FALSE]
    L2[, Col0] = 0
    L2[, Col1] = gX1pMX1 %*% X1pMX1
    L2[, Col2] = gX1pMX1 %*% X1pM %*% X2
    L2[abs(L2) < eps] = 0
    DF1 = NROW(L2[!apply(L2, 1, function(x) all(abs(x) < eps)), , drop=FALSE])

    R1 = t(M0[Col1, , drop=FALSE])
    R2 = t(M0[Col2, , drop=FALSE])
    L3 = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
    if (NROW(L3) > 0) {
      L3 = pivotJ(L3, Col0, clear=TRUE)
      L3 = L3[!apply(L3, 1, function(x) all(abs(x) < eps)), , drop=FALSE]
    }
    DF2 = NROW(L3)

    if (DF2 != DF1) {
      if (!is.null(Col2)) {
        Ms = pivotJ(M0[c(Col1, Col2), , drop=FALSE], Col0, clear=TRUE)
        Ms = pivotJ(Ms, Col1, clear=FALSE)

        fbazr = apply(abs(Ms[, Col1, drop=FALSE]), 1, max) < eps
        bazr = which(fbazr)
        bnazr = which(!fbazr)

        if (length(bnazr) == 0) {
          L3 = NULL
        } else if (length(bnazr) <= length(Col1)) {
          R1 = t(Ms[bnazr, , drop=FALSE])
          R2 = t(Ms[bazr, , drop=FALSE])
          L3 = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
          rownames(L3) = (rownames(M0)[Col1])[1:NROW(L3)]
        } else {
          R1 = t(M0[Col1, , drop=FALSE])
          R2 = t(M0[Col2, , drop=FALSE])
          L3 = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
        }
      } else {
        L3 = pivotJ(M0[Col1, , drop=FALSE], Col1, clear=FALSE)
        L3 = pivotJ(L3, Col0, clear=TRUE)
      }
      if (NROW(L3) > 0) {
        L3 = L3[!apply(L3, 1, function(x) all(abs(x) < eps)), , drop=FALSE]
      }
      DF3 = NROW(L3)
    } else {
      DF3 = DF1
    }

    if (DF3 == 0) L2 = matrix(0, nrow(L2), ncol(L2))
    re2 = rbind(re2, L2)
  }

  rownames(re2) = paste0("L", 1:NCOL(re2))
  re2[abs(re2) < eps2] = 0
  return(re2) # Do not use zapsmall !
}
