TTEST = function(x, y, conf.level=0.95)
{
  x = x[is.finite(x)]
  y = y[is.finite(y)]
  return(mtest(mean(x), sd(x), length(x), mean(y), sd(y), length(y), conf.level=conf.level))
}
