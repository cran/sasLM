G2SWEEP = function(A, Augmented=FALSE, eps=1e-8)
{
  p = nrow(A)
  p0 = ifelse(Augmented, p - 1, p)
  r = 0
  for (k in 1:p0) {
    d = A[k,k]
#    CSS = crossprod(A[k,] - mean(A[k,]))
#    DminK = ifelse(CSS > 0, eps*CSS, eps)
    if (abs(d) < eps) { A[k,] = 0 ; A[,k] = 0 ; next }
    A[k,] = A[k,]/d
    r = r + 1
    for (i in 1:p) {
      if (i != k) {
        c0 = A[i,k] ; 
        A[i,] = A[i,] - c0*A[k,] ; 
        A[i,k] = -c0/d
      }
    }
    A[k,k] = 1/d
  }
  attr(A, "rank") = r
  return(A)
}
