EMS = function(Formula, Data, Type=3, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  x = ModelMatrix(Formula, Data)
  Terms = labels(terms(x))
  nTerm = length(Terms)

  L0 = switch(as.character(Type),
              "1" = e1(Formula, Data, eps=eps),
              "2" = e2(Formula, Data, eps=eps),
              "3" = e3(Formula, Data, eps=eps),
              stop(paste("Type", Type, "is not supported!")))

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
      for (j in 1:nTerm) Res[i,j] = sum(diag(xC[x$assign==j, x$assign==j, drop=FALSE]))/NROW(L)
    } else {
#      Res[1,] = NA
      Res[i,] = 0
#      Res[i,i] = 1
    }
  }

  return(round(Res, abs(log10(eps)) + 1))
}
