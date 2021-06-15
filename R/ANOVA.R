ANOVA = function(Formula, Data, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")
  
  r1 = lfit(x, y, eps=eps)
  T1 = SS(x, r1, e1(Formula, Data, eps=eps))

  x2 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  r2 = lfit(x2, y, eps=eps)
  T2 = SS(x2, r2, e2(Formula, Data, eps=eps))[rownames(T1),,drop=FALSE]
  class(T2) = "anova"
  T3 = SS(x2, r2, e3(Formula, Data), eps=eps)[rownames(T1),,drop=FALSE]
  if (sum(T3[,"Df"], na.rm=TRUE) != sum(T1[,"Df"], na.rm=TRUE)) attr(T3, "heading") = "CAUTION: Singularity Exists !"
  class(T3) = "anova"

  if ("(Intercept)" %in% colnames(x$X)) {
    SST = crossprod(y - mean(y))
  } else {
    SST = crossprod(y)
  }

  ANOVA = sumANOVA(r1, T1=NULL, SST, nrow(x$X), rownames(attr(terms(x),"factors"))[1])

  Result = list(ANOVA=ANOVA, 'Type I'=T1, 'Type II'=T2, 'Type III'=T3)
  return(Result)
}
