GLM = function(Formula, Data, eps=1e-8)
{
  y = model.frame(Formula, Data)[,1]
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  r1 = lfit(x, y)
  T1 = SS(x, r1, e1(Formula, Data))

  x2 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  r2 = lfit(x2, y)
  T2 = SS(x2, r2, e2(Formula, Data))[rownames(T1),,drop=FALSE]
  class(T2) = "anova"
  T3 = SS(x2, r2, e3(Formula, Data))[rownames(T1),,drop=FALSE]
  if (sum(T3[,"Df"], na.rm=TRUE) != sum(T1[,"Df"], na.rm=TRUE)) attr(T3, "heading") = "CAUTION: Singularity Exists !"
  class(T3) = "anova"

  ANOVA = sumANOVA(r1, T1=NULL, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x),"factors"))[1])
  Parameter = sumREG(r1, x$X)

  Result = list(ANOVA=ANOVA, 'Type I'=T1, 'Type II'=T2, 'Type III'=T3, Parameter=Parameter)
  return(Result)
}
