SS = function(x, rx, L, eps=1e-8)
{
#  eps2 = min(abs(L[L != 0]))*eps
  eps2 = min(abs(rx$coefficients[rx$coefficients != 0]))*eps
  Labels = labels(terms(x))
  nLabel = length(Labels)
  ColNames = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  T1 = matrix(NA, nrow=nLabel, ncol=length(ColNames))
  dimnames(T1) = list(Labels, ColNames)

  for (i in 1:nLabel) {
    Li = L[x$termIndices[[Labels[i]]], , drop=FALSE]
    Li = Li[!apply(Li, 1, function(x) all(abs(x) < eps2)), , drop=FALSE]
    if (NROW(Li) > 0) {
      iiv = G2SWEEP(Li %*% rx$g2 %*% t(Li), Augmented=F, eps=eps2)
#      T1[i, "Df"] = NROW(Li)         # less accurate
      Lb.i = Li %*% rx$coefficients
      T1[i, "Df"] = attr(iiv, "rank") # more accurate
      T1[i, "Sum Sq"] = as.vector(t(Lb.i) %*% iiv %*% Lb.i)
    } else {
      T1[i, "Df"] = 0
      T1[i, "Sum Sq"] = NA
    }
  }

  T1[, "Mean Sq"] = T1[, "Sum Sq"]/T1[, "Df"]
  if (rx$DFr > 0) {
    T1[, "F value"] = T1[, "Mean Sq"]/(rx$SSE/rx$DFr)
    T1[, "Pr(>F)"] = 1 - pf(T1[, "F value"], T1[, "Df"], rx$DFr)
  } else {
    T1[, "F value"] = NA
    T1[, "Pr(>F)"] = NA
  }
  class(T1) = "anova"

  return(T1)
}
