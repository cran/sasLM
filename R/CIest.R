CIest = function(Formula, Data, Term, Contrast=c(-1, 1), Alpha=0.10) 
{
  x = ModelMatrix(Formula, Data)
  if (!(Term %in% attr(x$terms, "term.labels"))) stop(paste("There is no term: ", Term))

  colIndex = x$termIndices[Term][[1]]
  if (length(Contrast) != length(colIndex)) stop("Constrast length is not appropriate!")

  y = model.frame(Formula, Data)[,1]
  r2 = lfit(x, y)

  L = rep(0, length(r2$coefficients))
  L[colIndex] = Contrast

  Diff = est(t(L), r2)
  t0 = qt(1 - Alpha/2, r2$DFr)
  PE = Diff[1, "Estimate"]
  CI = PE + c(-1, 1)*t0*Diff[1, "Std. Error"]

  Result = c(PE, CI)
  names(Result) = c("Point Estimate", "Lower Limit", "Upper Limit")

  return(Result)
}
