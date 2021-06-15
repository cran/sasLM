SLICE = function(Formula, Data, mTerm, sTerm)
{
  if (!attr(terms(Formula, data=Data), "response")) stop("Dependent variable should be provided!")
  x = ModelMatrix(Formula, Data)

  lLabel = attr(terms(x), "term.labels")
  if (!(mTerm %in% lLabel)) stop("mTerm(mean term) is not found in the formula!")
  if (!(sTerm %in% lLabel)) stop("sTerm(slice term) is not found in the formula!")
  
  r1 = aov3(Formula, Data)
  MSE = r1["RESIDUALS", "Mean Sq"]
  DFr = r1["RESIDUALS", "Df"]   
  
  Levels = unique(Data[,sTerm])
  nLevel = length(Levels)
  Res = NULL
  for (i in 1:nLevel) {
    r1s = aov3(Formula, Data[Data[,sTerm] == Levels[i],])[paste0(" ", mTerm),]
    Res = rbind(Res, r1s)
  }
  Res[,"F value"] = Res[,"Mean Sq"]/MSE
  Res[,"Pr(>F)"] = 1 - pf(Res[,"F value"], Res[,"Df"], DFr)
  rownames(Res) = Levels
  class(Res) = c("anova", class(Res))
  return(Res)
}
