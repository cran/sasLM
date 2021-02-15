LCL = function(x, conf.level=0.95)
{
  if (!(is.numeric(x) & is.vector(x))) stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  mean(x) + qt(0.5 - conf.level/2, n - 1)*sd(x)/sqrt(n)
}
