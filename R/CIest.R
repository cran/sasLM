CIest = function(Formula, Data, Term, Contrast, conf.level=0.95) 
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  if (missing(Term) & length(labels(terms(x))) == 1) Term = labels(terms(x))
  if (!(Term %in% attr(x$terms, "term.labels"))) stop(paste("There is no term: ", Term))

  colIndex = x$termIndices[Term][[1]]
  if (length(Contrast) != length(colIndex)) stop("Contrast length is not appropriate!")

  r2 = lfit(x, y)

  L = rep(0, length(r2$coefficients))
  L[colIndex] = Contrast

  Res = est(t(L), x$X, r2, conf.level=conf.level)
#  printCoefmat(Res)
#  invisible(Res)
  class(Res) = "anova"
  return(Res)
}
