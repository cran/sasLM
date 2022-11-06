UCL = function(y, conf.level=0.95) 
{
  y = as.numeric(y)  
  if (!(is.numeric(y) & is.vector(y))) return(NA) # stop("Input should be a numeric vector!")
  y = y[!is.na(y)]
  n = length(y)
  if (n < 2) return(NA) # stop("The length of the input vector should be larger than 1.")
  mean(y) + qt(0.5 + conf.level/2, n - 1)*sd(y)/sqrt(n)
}

