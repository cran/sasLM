GLM = function(Formula, Data, BETA=FALSE, EMEAN=FALSE, conf.level=0.95, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data, KeepOrder=TRUE)
  mf0 = model.frame(Formula, Data)
  y = mf0[, 1]
  yName = names(mf0)[1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  r1 = lfit(x, y, eps=eps)
  T1 = SS(x, r1, e1(Formula, Data, eps=eps), eps=eps)

  x2 = ModelMatrix(Formula, Data, KeepOrder=FALSE)
  r2 = lfit(x2, y, eps=eps)
  T2 = SS(x2, r2, e2(Formula, Data, eps=eps), eps=eps)[rownames(T1),,drop=FALSE]
  class(T2) = "anova"
  T3 = SS(x2, r2, e3(Formula, Data, eps=eps), eps=eps)[rownames(T1),,drop=FALSE]
  if (sum(T3[,"Df"], na.rm=TRUE) != sum(T1[,"Df"], na.rm=TRUE)) attr(T3, "heading") = "CAUTION: Singularity Exists !"
  class(T3) = "anova"

  if ("(Intercept)" %in% colnames(x$X)) {
    SST = crossprod(y - mean(y))
  } else {
    SST = crossprod(y)
  }

  ANOVA = sumANOVA(r1, T1=NULL, SST, nrow(x$X), rownames(attr(terms(x),"factors"))[1])

  Rsq = ANOVA["MODEL", "Sum Sq"]/ANOVA[3, "Sum Sq"] # Corrected total or uncorrected total
  RMSE = sqrt(ANOVA["RESIDUALS", "Mean Sq"])
  MeanY = mean(y, na.rm=T)
  CV = 100*RMSE/MeanY
  Fit = data.frame(Rsq, CV, RMSE, MeanY)
  colnames(Fit) = c("R-square", "Coef Var", "Root MSE", paste(yName, "Mean"))
  rownames(Fit) = ""

  Result = list(ANOVA=ANOVA, Fitness=Fit, 'Type I'=T1, 'Type II'=T2, 'Type III'=T3)
  iNext = 6
  
  if (BETA) {
    Result[[iNext]] = sumREG(r1, x$X)
    names(Result)[iNext] = "Parameter"
    iNext = iNext + 1
  }

  if (EMEAN) {
    Result[[iNext]] = lsm0(x2, r2, Formula, Data, conf.level = conf.level)
    names(Result)[iNext] = "Expected Mean"
  }

  return(Result)
}
