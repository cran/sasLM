Cor.test = function(Data, conf.level=0.95)
{
  nCol = ncol(Data)
  vNum = rep(FALSE, nCol)
  for (i in 1:nCol) vNum[i] = is.numeric(Data[,i])
  Data = Data[,vNum]

  nc0 = ncol(Data)
  if (nc0 < 2) stop("Input should have more than 1 column!")
  for (i in 1:nc0) if (!is.numeric(Data[,i])) stop("This is only for numeric columns!")
  
  cn0 = colnames(Data)
  mC = cor(Data, use="pairwise.complete.obs")

  m2 = matrix(as.numeric(!is.na(Data)), ncol=nc0)
  colnames(m2) = cn0
  x = model.matrix(~ . - 1, data.frame(m2))
  n = crossprod(x)

  Df = n - 2
  Tval = mC*sqrt(Df/(1 - mC^2))
  Pval = 2*pt(-abs(Tval), Df)
  Zs = atanh(mC)
  DL = qnorm(0.5 + conf.level/2)/sqrt(n - 3)
  LL = tanh(Zs - DL)
  UL = tanh(Zs + DL)
  
  ColNames = c("Estimate", "Lower CL", "Upper CL", "t value", "Df", "Pr(>|t|)")
  nr =  nc0*(nc0 - 1)/2
  nc = length(ColNames)
  Res = matrix(nrow=nr, ncol=nc)
  colnames(Res) = ColNames
  rownames(Res) = 1:nr
  cr = 1
  for (i in 1:nc0) {
    for (j in (i + 1):nc0) {
      Res[cr, 1] = mC[i, j]
      Res[cr, 2] = LL[i, j]
      Res[cr, 3] = UL[i, j]
      Res[cr, 4] = Tval[i, j]
      Res[cr, 5] = Df[i, j]
      Res[cr, 6] = Pval[i, j]
      rownames(Res)[cr] = paste0(cn0[i], ":", cn0[j])
      cr = cr + 1
    }
    if (cr > nr) break
  }
  printCoefmat(Res)
  invisible(Res)
}
