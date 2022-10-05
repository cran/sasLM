pResD = function(rx, Title=NULL)
{
  X = model.matrix(rx)
  n = nrow(X)
  p = qr(X)$rank
  yhat = rx$fitted
  e = rx$residuals
  y = yhat + e
  h = hatvalues(rx)
  Rst = rstandard(rx)
  CooksD = cooks.distance(rx)

  oPar = par(mfrow=c(3, 3))
  if (!is.null(Title)) par(oma=c(1, 1, 2, 1))

  rngY0 = c(min(e, -3), max(e, 3))
  plot(e ~ yhat, ylim=rngY0, xlab="Predicted Value", ylab="Residual", main="PLOT A", cex.main=0.9)
  abline(h=0, lty=3)

  rngY1 = c(min(Rst, -2), max(Rst, 2))
  plot(Rst ~ yhat, ylim=rngY1, xlab="Predicted Value", ylab="RStudent", main="PLOT B", cex.main=0.9)
  abline(h=c(-2, 2), lty=3)

  plot(Rst ~ h, ylim=rngY1, xlab="Leverage", ylab="RStudent", main="PLOT C", cex.main=0.9)
  abline(h=c(-2, 2), v=2*p/n, lty=3)

  qqnorm(e, main="PLOT D", cex.main=0.9)
  qqline(e, lty=3)

  plot(y ~ yhat, xlab="Predicted Value", ylab="Observed Y value", main="PLOT E", cex.main=0.9)
  abline(a=0, b=1, lty=3)

  plot(CooksD, xlab="Observation", ylab="Cook's D", type="h", main="PLOT F", cex.main=0.9)
  points(CooksD)
  abline(h=c(0, 4/n), lty=3)

  hist(e, freq=F, breaks=sqrt(n), xlab="Residual", ylab="Density", main="PLOT G", cex.main=0.9)
  lines(density(e))

  rngY = range(yhat - mean(y))
  qqnorm(yhat - mean(y), ylim=rngY, main="Plot H", cex.main=0.9)
  qqline(yhat - mean(y), lty=3)

  qqnorm(e, ylim=rngY, xlab="Residual", main="Plot I", cex.main=0.9)

  if (!is.null(Title)) title(Title, outer=TRUE)
  par(oPar)
}
