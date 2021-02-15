lr0 = function(Formula, Data)
{
  mf = model.frame(Formula, Data)
  cn = colnames(mf)[-1]
  ni = ncol(mf) - 1
  my = mean(mf[,1])
  SST = as.numeric(crossprod(mf[,1] - my))

  Res = matrix(nrow=ni, ncol=6)
  colnames(Res) = c("Intercept", "SE(Intercept)", "Slope", "SE(Slope)", "Rsq", "Pr(>F)")
  rownames(Res) = cn

  for (i in 1:ni) {
    tr = summary(lm(mf[,1] ~ mf[,i + 1]))
    Res[i, 1:2] = tr$coefficients[1, 1:2]
    Res[i, 3:4] = tr$coefficients[2, 1:2]
    Res[i, 5] = tr$r.squared
    Res[i, 6] = tr$coefficients[2, 4]
  }
  printCoefmat(Res)
  invisible(Res)
}
