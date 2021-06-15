pB = function(Formula, Data, Resol=300, conf.level=0.95, lx, ly, ...)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  
  mf = model.frame(Formula, Data)
  if (!is.numeric(mf[,1])) stop("Dependent variable should be numeric!")

  nc = ncol(mf)
  if (nc != 2) stop("Only one dependent and one independent variable is supported!")
  
  r1 = lm(Formula, Data)
  
  if ("xlim" %in% names(list(...))) {
    rangeX = list(...)$xlim
  } else {
    rangeX = range(mf[,2])
  }

  x1 = seq(rangeX[1], rangeX[2], length.out=(Resol + 1))
  d1 = data.frame(x1)
  colnames(d1) = colnames(mf)[2]
  p1 = predict(r1, d1, interval="prediction", level=conf.level)
  c1 = predict(r1, d1, interval="confidence", level=conf.level)

  if ("ylim" %in% names(list(...))) {
    plot(Formula, Data, ...)
  } else {  
    rangeY = range(p1)
    plot(Formula, Data, ylim=rangeY, ...)
  }

  matlines(x1, c1, lty=c(1,2,2), col=c(1,2,2))
  matlines(x1, p1, lty=c(1,4,4), col=c(1,4,4))
  if ((missing(lx) | missing(ly)) & coef(r1)[2] > 0) {
    lx = rangeX[1] + 0.7*(rangeX[2] - rangeX[1])
    ly = rangeY[1] + 0.1*(rangeY[2] - rangeY[1])
  } else {
    lx = rangeX[1] + 0.7*(rangeX[2] - rangeX[1])
    ly = rangeY[2]    
  }
  legend(lx, ly, c("Confidence Interval", "Prediction Interval"), lty=c(2,4), col=c(2,4), cex=0.75)
}
