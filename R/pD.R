pD = function(rx, Title=NULL)
{
  yhat = fitted(rx)
  y = yhat + residuals(rx)
  sRes = rstandard(rx)

#  dev.new()
  oPar = par(mfrow=c(2,2))
  if(!is.null(Title)) par(oma=c(1, 1, 2, 1))

  y0 = min(c(y, yhat))
  y1 = max(c(y, yhat))
  plot(yhat, y, xlim=c(y0, y1), ylim=c(y0, y1), xlab=expression(hat("y")), ylab="y", pch=16, main="")
  abline(a=0, b=1, lty=3)
  mtext("Observed vs. Fitted Values", line=0.5)

  rM = max(abs(sRes))
  plot(yhat, sRes, ylim=c(-rM, rM), xlab=expression(hat("y")), ylab="Standardized Residual", pch=16, main="")
  lines(lowess(sRes ~ yhat), col="red")
  abline(h=0, lty=3)
  mtext("Standardized Residual vs. Fitted", line=0.5)

  yM = max(c(hist(sRes, plot=FALSE)$density, density(sRes)$y))
  hist(sRes, freq=FALSE, ylim=c(0, yM), xlab="Standardized Residual", ylab="Density", main="")
  lines(density(sRes))
  mtext("Distribution of Standardized Res.", line=0.5)

  qqnorm(sRes, ylab="Standardized Residual", main="")
  qqline(sRes, lty=3)
  mtext("Q-Q Plot of Standardized Res.", line=0.5)

  if(!is.null(Title)) title(Title, outer=TRUE)
  par(oPar)
}

