G2SWEEP = function(A, Augmented=FALSE, eps=1e-8)
{
#  eps2 = 1e-14
  idx = abs(diag(A)) > eps
  p = sum(idx, na.rm=T)
  p0 = ifelse(Augmented, p - 1, p)
  if (p == 0 | p0 < 1) { A[, ] = 0 ; attr(A, "rank") = 0 ; return(A) }
  B = A[idx, idx, drop=F]

#  if (!Augmented) {
#    s = apply(B, 1, function(x) (max(abs(x))*min(abs(x[x!=0])))^(1/4))
#    B = as.matrix(apply(B, 1, function(x) x/s))
#    B = as.matrix(apply(B, 2, function(x) x/s))
#  }

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

#  if (!Augmented) {
#    B = as.matrix(apply(B, 1, function(x) x/s))
#    B = as.matrix(apply(B, 2, function(x) x/s))
#    B[abs(B) < eps2] = 0
#  }

  A[!idx, !idx] = 0
  A[idx, idx] = B
  attr(A, "rank") = r

  return(A)
}

