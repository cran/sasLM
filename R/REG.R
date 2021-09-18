REG = function(Formula, Data, eps=1e-8, summarize=TRUE)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  if ("Complete" %in% names(CheckAlias(Formula, Data))) {
    warning("Complete aliased variable(s) exist(s)!")
    eps = 1e-5
  }

  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  if (summarize) {
    return(sumREG(lfit(x, y, eps=eps), x$X))
  } else {
    return(lfit(x, y, eps=eps))
  }
}
