e1 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, NOINT=FALSE, KeepOrder=TRUE)
  A = crossprod(x$X)

  nr = nrow(A)
  for (i in 1:(nr - 1)) {
    if (abs(A[i,i]) < eps) { A[i,] = 0 ; A[i:nr,i] = 0 ; next }        
    for (j in (i + 1):nr)  A[j,] = A[j,] - A[j,i]/A[i,i]*A[i,]
    if (abs(A[i,i]) > eps) A[i,] = A[i,]/A[i,i]
  }
  colnames(A) = paste0("L", 1:nrow(A))
  return(round(A, 11))
}
