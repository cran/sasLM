is.cor = function(m, eps=1e-16)
{
  nc = ncol(m)
  if (nc != nrow(m)) stop("Input should be a sqaure matrix!")
  for (i in 1:nc) if (!is.numeric(m[,i])) stop("All columns should be numeric!")
  for (i in 1:nc) {
    for (j in i:nc) {
      if (j == i) {
        if (abs(m[i, i] - 1) > eps) return(FALSE)
      } else {
        if (abs(m[i, j] > 1)) return(FALSE)
        if (abs(m[i, j] - m[j, i]) > eps) return(FALSE)
      }
    }
  }
  return(TRUE)
}
