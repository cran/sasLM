T3MS = function(Formula, Data, L0, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  y = model.frame(Formula, Data)[, 1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x = ModelMatrix(Formula, Data)
  Terms = labels(terms(x))
  nTerm = length(Terms)

  if (missing(L0)) L0 = e3(x, eps=eps)

  r0 = lfit(x, y)
  b = r0$coefficients

  Res = matrix(nrow=nTerm, ncol=nTerm)
  rownames(Res) = Terms
  colnames(Res) = Terms
  for (i in 1:nTerm) {
    L = L0[(x$assign == i) & (abs(b) > eps), , drop=FALSE]
    if (NROW(L) > 0) {
      xC = t(L) %*% G2SWEEP(L %*% G2SWEEP(crossprod(x$X)) %*% t(L)) %*% L
#      M = qr.solve(t(chol(L %*% G2SWEEP(crossprod(x$X)) %*% t(L)))) # Frequent crash
#      xC = crossprod(M %*% L)
      for (j in 1:nTerm) Res[i,j] = sum(diag(xC[x$assign==j, x$assign==j]))/NROW(L)
    } else {
      Res[i, ] = 0
    }
  }

  Res[abs(Res) < eps] = 0
  return(Res)
}
