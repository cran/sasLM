WhiteTest = function(rx)
{
  if (!("lm" %in% class(rx))) stop("Input should be a result of lm")
  if (is.null(weights(rx))) { w = rep(1, nobs(rx))
  } else { w = weights(rx) }
  X = sqrt(w)*model.matrix(rx)
  we2 = w*resid(rx)^2     # we2 = weight * e^2
  return(WhiteH(X, we2))
}
