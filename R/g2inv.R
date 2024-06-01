g2inv = function(A, eps=1e-8)
{
  idx = abs(diag(A)) > eps
  p = sum(idx, na.rm=T)
  M = matrix(0, nrow=NCOL(A), ncol=NROW(A))
  if (p == 0) {attr(M, "rank") = 0 ; return(M) }
  B = A[idx, idx, drop=F]
  r = 0
  for (k in 1:p) {
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
  M[1:r, 1:r] = B[1:r, 1:r]
  attr(M, "rank") = r
  return(M)
}
