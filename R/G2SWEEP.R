G2SWEEP = function(A, Augmented=FALSE, eps=1e-8)
{
  idx = abs(diag(A)) > eps
  p = sum(idx)
  if (p == 0) {A[,] = 0 ; return(A)}
  B = A[idx, idx, drop=F]

  p0 = ifelse(Augmented, p - 1, p)
  r = 0
  for (k in 1:p0) {
    d = B[k, k]
    if (abs(d) < eps) { B[k, ] = 0 ; B[, k] = 0 ; next }
    B[k, ] = B[k, ]/d
    r = r + 1
    for (i in 1:p) {
      if (i != k) {
        c0 = B[i, k]
        B[i, ] = B[i, ] - c0*B[k, ]
        B[i, k] = -c0/d
      }
    }
    B[k, k] = 1/d
  }

  A[!idx, !idx] = 0
  A[idx, idx] = B
  attr(A, "rank") = r
  return(A)
}
