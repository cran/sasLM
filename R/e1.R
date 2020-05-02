e1 = function(Formula, Data, eps=1e-8)
{
  x = ModelMatrix(Formula, Data, NOINT=FALSE, KeepOrder=TRUE)
  XpX = crossprod(x$X)

  nr = nrow(XpX)
  for (i in 1:(nr - 1)) {
    if (abs(XpX[i,i]) < eps) { XpX[i,] = 0 ; XpX[i:nr,i] = 0 ; next }        
    for (j in (i + 1):nr)  XpX[j,] = XpX[j,] - XpX[j,i]/XpX[i,i]*XpX[i,]
    if (abs(XpX[i,i]) > eps) XpX[i,] = XpX[i,]/XpX[i,i]
  }
  rownames(XpX) = paste0("L", 1:nrow(XpX))
  return(XpX)
}
