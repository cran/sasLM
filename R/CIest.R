CIest = function(Formula, Data, Term, Contrast, conf.level=0.95) 
{
  x = ModelMatrix(Formula, Data)
  if (!(Term %in% attr(x$terms, "term.labels"))) stop(paste("There is no term: ", Term))

  colIndex = x$termIndices[Term][[1]]
  if (length(Contrast) != length(colIndex)) stop("Constrast length is not appropriate!")

  y = model.frame(Formula, Data)[,1]
  r2 = lfit(x, y)

  L = rep(0, length(r2$coefficients))
  L[colIndex] = Contrast

  Res = est(t(L), x$X, r2, conf.level=conf.level)
  printCoefmat(Res)
  invisible(Res)
}
