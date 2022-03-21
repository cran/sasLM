UCL = function(x, conf.level=0.95) 
{
  x = as.numeric(x)  
  if (!(is.numeric(x) & is.vector(x))) return(NA) # stop("Input should be a numeric vector!")
  x = x[!is.na(x)]
  n = length(x)
  if (n < 2) return(NA) # stop("The length of the input vector should be larger than 1.")  
  mean(x) + qt(0.5 + conf.level/2, n - 1)*sd(x)/sqrt(n)
}

