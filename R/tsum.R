tsum = function(Formula=NULL, Data=NULL, ColNames=NULL, MaxLevel=30, ...)
{
  if (is.data.frame(Formula) | is.matrix(Formula)) {
    Data = as.data.frame(Formula)
    Formula = NULL
  } else if (is(Formula, "vector") & is.numeric(Formula)) {
    Data = as.data.frame(Formula)
    colnames(Data) = as.character(as.list(match.call())$Formula)
    Formula = NULL
  } else if (is(Formula, "formula") & !is.data.frame(Data)) {
    if (is.null(Data)) { stop("Data should be provided with Formula.")
    } else { Data = as.data.frame(Data) }
  }

  if (is(Formula, "formula")) {
    mf = model.frame(Formula, Data)
    nc = ncol(mf)
    if (nc > 4) stop("Too many independent variables.")
    cn = colnames(mf)
    yn = cn[1]
    if (nc > 1) {
      cn = cn[-1]
      nx = nc - 1
      vc = vector(length=nx)
      for (i in 1:nx) {
        vc[i] = length(unique(mf[, i + 1]))
        if (vc[i] > MaxLevel) stop("Too many levels. Increase MaxLevel if you want!")
      }
      cn = cn[order(vc)]
      Res = switch(as.character(nx),
                 '1' = tsum1(mf, yn, cn[1], ...),
                 '2' = tsum2(mf, yn, cn[1], cn[2], ...),
                 '3' = tsum3(mf, yn, cn[1:2], cn[3], ...))
    } else {
      Res = tsum0(mf, yn, ...)
    }
  } else {
    if (is.null(ColNames)) ColNames = colnames(Data)
    nCol = length(ColNames)
    Res = list()
    cr = 1
    for (i in 1:nCol) {
      if (is.numeric(Data[,ColNames[i]])) {
        Res[[cr]] = tsum0(Data, ColNames[i], ...)
        names(Res)[cr] = ColNames[i]
        cr = cr + 1
      }
    }
    Res = as.data.frame(Res)
  }
  return(Res)
}
