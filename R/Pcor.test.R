Pcor.test = function(Data, x, y)
{
  d = Data[,c(x, y)]
  d = d[complete.cases(d),]
  cn0 = colnames(d)
  n = nrow(d)
  Df = n - length(y) - 2
  m = cor(d)
  mC = cov2cor(m[x, x] - m[x, y] %*% solve(m[y, y], t(m[x, y])))
  Tval = mC * sqrt(Df/(1 - mC^2))
  Pval = 2*pt(-abs(Tval), Df)
  
  nc0 = length(x)
  nr =  nc0*(nc0 - 1)/2
  ColNames = c("Estimate", "Df", "t value", "Pr(>|t|)")
  nc = length(ColNames)
  Res = matrix(nrow=nr, ncol=nc)
  colnames(Res) = ColNames
  rownames(Res) = 1:nr
  cr = 1
  for (i in 1:nc0) {
    for (j in (i + 1):nc0) {
      Res[cr, 1] = mC[i, j]
      Res[cr, 2] = Df
      Res[cr, 3] = Tval[i, j]
      Res[cr, 4] = Pval[i, j]
      rownames(Res)[cr] = paste0(cn0[i], ":", cn0[j])
      cr = cr + 1
    }
    if (cr > nr) break
  }
  printCoefmat(Res)
  invisible(Res)
}
