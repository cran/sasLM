SLICE = function(Formula, Data, Term, By)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }
  x = ModelMatrix(Formula, Data)
  lLabel = attr(terms(x), "term.labels")
  if (!(Term %in% lLabel)) stop("Term(mean term) is not found in the formula!")
  if (!(By %in% lLabel)) stop("By(slice term) is not found in the formula!")
  
  r1 = aov3(Formula, Data)
  MSE = r1["RESIDUALS", "Mean Sq"]
  DFr = r1["RESIDUALS", "Df"]   
  
  Levels = unique(Data[, By])
  nLevel = length(Levels)
  Res = NULL
  for (i in 1:nLevel) {
    fNew = formula(paste0(terms(Formula)[[2]], " ~ ", Term))
    dNew = Data[Data[, By] == Levels[i], ]
    r1s = aov3(fNew, dNew)[paste0(" ", Term), ]
    Res = rbind(Res, r1s)
  }
  Res[, "F value"] = Res[, "Mean Sq"]/MSE
  Res[, "Pr(>F)"] = 1 - pf(Res[, "F value"], Res[, "Df"], DFr)
  rownames(Res) = Levels
  class(Res) = c("anova", class(Res))
  return(Res)
}
