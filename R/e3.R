e3 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, KeepOrder=FALSE)
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

    Col1 = x$termIndices[Label1][[1]] # direct effect
    Col2 = NULL                       # containing effect
    if (length(Label2) > 0) {
      for (j in 1:length(Label2)) Col2 = c(Col2, x$termIndices[Label2[j]][[1]])
    }
    Col0 = setdiff(1:nc, c(Col1, Col2)) # unrelated effect

    if (!is.null(Col2)) {
## Get Degree of Freedom (DF) using e2 method      
      X0 = x$X[,Col0]
      X1 = x$X[,Col1]
      X2 = x$X[,Col2]
      Mx = X0 %*% G2SWEEP(crossprod(X0)) %*% t(X0)
      M0 = diag(NCOL(Mx)) - Mx

      X1pM  = crossprod(X1, M0)
      X1pMX1 = X1pM %*% X1
      gX1pMX1 = G2SWEEP(X1pMX1)

      Lx = M[x$termIndices[Label1][[1]], , drop=FALSE]
      Lx[,Col0] = 0
      Lx[,Col1] = gX1pMX1 %*% X1pMX1
      Lx[,Col2] = gX1pMX1 %*% X1pM %*% X2
      Lx = Lx[!apply(Lx, 1, function(x) all(abs(x) < eps)), , drop=FALSE]
      DF1 = NROW(Lx)   
## End of getting DF

## Try conventional e3 method first
      R1 = t(M[Col1, , drop=FALSE])
      R2 = t(M[Col2, , drop=FALSE])
      tL = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
      if (NROW(tL) > 0) {
        tL = pivotJ(tL, Col0, clear=TRUE)
        tL = tL[!apply(tL, 1, function(x) all(abs(x) < eps)), , drop=FALSE]
      }
      DF2 = NROW(tL)   
## End of getting conventional e3

## If the conventional e3 does not work, i.e. rare case (2 cases / 194 cases)
      if (DF2 != DF1) {
        Ms = pivotJ(M[c(Col1, Col2), , drop=FALSE], Col0, clear=TRUE)
        Ms = pivotJ(Ms, Col1, clear=FALSE)

        fbazr = apply(abs(Ms[, Col1, drop=FALSE]), 1, max) < eps
        bazr = which(fbazr)
        bnazr = which(!fbazr)

        if (length(bnazr) == 0) {
          tL = NULL
        } else if (length(bnazr) <= length(Col1)) {
          R1 = t(Ms[bnazr, , drop=FALSE])
          R2 = t(Ms[bazr, , drop=FALSE])
          tL = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
          rownames(tL) = (rownames(M)[Col1])[1:NROW(tL)]
        } else {
          R1 = t(M[Col1, , drop=FALSE])
          R2 = t(M[Col2, , drop=FALSE])
          tL = t(R1 - R2 %*% G2SWEEP(crossprod(R2)) %*% crossprod(R2, R1))
        }
      }
## End of rare case

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
