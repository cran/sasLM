lfit = function(x, y, eps=1e-8)
{
  if (NCOL(y) > 1) stop("This support only 1 column of y.")
  nc = ncol(x$X)
  XpY = crossprod(x$X, y)
  aXpX = rbind(cbind(crossprod(x$X), XpY), cbind(t(XpY), crossprod(y)))
  ag2 = G2SWEEP(aXpX, Augmented=TRUE, eps=eps)

  b = ag2[1:nc,(nc + 1)]
  iXpX = ag2[1:nc, 1:nc]
  SSE = ag2[(nc + 1), (nc + 1)]
  DFr = nrow(x$X) - attr(ag2, "rank")
  
  Result = list(coefficients=b, g2=iXpX, rank=attr(ag2, "rank"), DFr=DFr, SSE=SSE)
  return(Result)
}

