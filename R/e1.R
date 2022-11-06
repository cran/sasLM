e1 = function(XpX, eps=1e-8)
{ # Forward Doolittle Method
  eps2 = eps/max(abs(XpX))
  nr = nrow(XpX)
  if (nr > 1) {
    for (i in 1:(nr - 1)) {
      if (abs(XpX[i, i]) < eps) { XpX[i, ] = 0 ; XpX[i:nr, i] = 0 ; next }
      for (j in (i + 1):nr)  XpX[j, ] = XpX[j, ] - XpX[j, i]/XpX[i, i]*XpX[i, ]
      if (abs(XpX[i, i]) > eps) XpX[i, ] = XpX[i, ]/XpX[i, i]
    }
    XpX[abs(XpX) < eps2] = 0
  } else {
    XpX[1, 1] = 1
  }
  rownames(XpX) = paste0("L", 1:nrow(XpX))
  return(XpX) # Do not use zapsmall !
}
