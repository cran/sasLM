lfit = function(x, y, eps=1e-8)
{
  if (NCOL(y) > 1) stop("This supports only y with 1 column.")
  nc = ncol(x$X)
  XpY = crossprod(x$X, y)
  XpX = crossprod(x$X)
  aXpX = rbind(cbind(XpX, XpY), cbind(t(XpY), crossprod(y)))
  ag2 = G2SWEEP(aXpX, Augmented=TRUE, eps=eps)

  b = ag2[1:nc,(nc + 1)]
  iXpX = ag2[1:nc, 1:nc, drop=FALSE]
  SSE = max(0, ag2[(nc + 1), (nc + 1)]) # SSE cannot be negative
  DFr = nrow(x$X) - attr(ag2, "rank")
  DFr2 = rep(DFr, nc)
  DFr2[diag(XpX) == 0] = NA

#  SST = as.numeric(crossprod(y - mean(y)))
#  R2 = 1 - SSE/SST
#  n = length(y)
#  R2ADJ = 1 - (1 - R2)*(n - 1)/DFr

  Result = list(coefficients=b, g2=iXpX, rank=attr(ag2, "rank"), DFr=DFr, SSE=SSE, DFr2=DFr2)
  return(Result)
}
