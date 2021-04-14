estmb = function(L, X, g2, eps=1e-8) # Estimability check
{ # Reference: Kennedy & Gentle. Statistical Computing (1980) p361 
  nc = ncol(L)
  if (nc != ncol(X) | nc != ncol(g2)) stop ("Matrix dimension mismatch!")
  L2 = L %*% g2 %*% crossprod(X)
  nr = nrow(L)
  vL = vector(length=nr)
  for (i in 1:nr) vL[i] = all(abs(L[i,] - L2[i,]) < eps)
  return(vL)
}

