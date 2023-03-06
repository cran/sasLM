RanTest = function(Formula, Data, Random, Type=3, eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  if (missing(Random)) {
    stop("At least one random factor should be provided!")
  }

  if (Type == 1) {
    x = ModelMatrix(Formula, Data, KeepOrder = TRUE)
  } else {
    x = ModelMatrix(Formula, Data, KeepOrder = FALSE)
  }
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  L0 = switch(as.character(Type),
              "1" = e1(crossprod(x$X), eps=eps),
              "2" = e2(x, eps=eps),
              "3" = e3(x, eps=eps),
              stop(paste("Type", Type, "is not supported!")))

  r1 = lfit(x, y, eps=eps)
  T1 = SS(x, r1, L0)
  T2 = rbind(T1, c(Df=r1$DFr, r1$SSE, r1$SSE/r1$DFr, NA, NA))
  rownames(T2) = c(rownames(T1), "MSE")

  RowNamesT1 = rownames(T1)
  strTerm = strsplit(RowNamesT1, ":")
  nError = 0
  Error = NULL
  for (i in 1:length(RowNamesT1)) {
    if (length(strTerm[[i]]) > 1 & length(intersect(Random, strTerm[[i]])) > 0) {
      nError = nError + 1
      Error = c(Error, RowNamesT1[i])
    }
  }

  T0 = EMS(Formula, Data, Type=Type, eps=eps)
  ColNames = strsplit(colnames(T0), ":")
  RowNames = strsplit(rownames(T0), ":")
  nCol = ncol(T0)

  W0 = matrix(0, ncol=nCol, nrow=nCol)
  dimnames(W0) = dimnames(T0)

  T5 = T0
  for (j in 1:nCol) {
    if (length(intersect(Random, ColNames[[j]])) == 0) T5[-j, j] = 0
  }

  for (i in 1:(nCol - 1)) {
    for (j in (i + 1):nCol) {
      if (T5[i, j] != 0) {
        W0[i, j] = T5[i, j]/T5[j, j]
        if (j < nCol) {
          for (k in (j + 1):nCol) {
            if (T5[j, k] != 0) T5[i, k] = T5[i, k] - W0[i, j]*T5[j, k]
          }
        }
      }
    }
  }
  W0 = cbind(W0, MSE = round(1 - apply(W0, 1, sum), 8))

  T4 = as.data.frame(round(T0, 4))
  for (j in 1:nCol) {
    if (length(intersect(Random, ColNames[[j]])) == 0) T4[T0[, j] > 0, j] = "Q"
  }
  T4[T4 == 0 | T4 == "0"] = ""
  fQ = ifelse(any(T4[upper.tri(T4)] == "Q"), TRUE, FALSE)
  T4 = cbind(T4, eps=1)

  ToTest = NULL
  for (i in 1:length(Error)) {
    ToTest = union(ToTest, setdiff(rownames(T0[which(T0[, Error[i]] > 0), , drop=F]), Error[i]))
  }
  nTest = length(ToTest)

  Result = list()
  if (nTest > 0) {
    for (i in 1:nTest) {
      v1 = T2[ToTest[i], "Mean Sq"]
      df1 = T2[ToTest[i], "Df"]

      df2 = T2[union(colnames(W0)[abs(W0[ToTest[i], ]) > eps], "MSE"), "Df"]
      v2 = T2[names(df2), "Mean Sq"]
      w2 = W0[ToTest[i], ][names(df2)]
      v4 = satt(vars=v2, dfs=df2, ws=w2)

      DF = c(df1, v4$Df)
      SS = c(T1[ToTest[i], "Sum Sq"], v4$Variance * v4$Df)
      MS = c(v1, v4$Variance)
      Fval = c(v1/v4$Variance, NA)
      Pval = c(1 - pf(Fval[1], df1, v4$Df), NA)

      T3 = cbind(DF, SS, MS, Fval, Pval)
      colnames(T3) = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
      rownames(T3) = c(ToTest[i], "Error")

      nErrLen = length(w2)
      wName = names(w2)
      strErr = paste0(round(w2[1], 4), "*", wName[1])
      for (j in 2:nErrLen) {
        strErr = paste0(strErr, ifelse(w2[[j]] >= 0, " + ", " - "), round(abs(w2[[j]]), 4), "*", wName[j])
      }
      attr(T3, "heading") = paste("Error:", strErr)
      class(T3) = "anova"
      Result[[i]] = T3
    }

    T2b = T2[setdiff(rownames(T2), ToTest), ]
    rownames(T2b) = c(rownames(T2b)[-nrow(T2b)], "RESIDUALS")
    class(T2b) = "anova"
    Result[[nTest + 1]] = T2b
    names(Result) = c(ToTest, "The Rest Terms")

    nNext = length(Result) + 1
    Result[[nNext]] = T4
    names(Result)[nNext] = "EMS"
  } else {
    rownames(T2)[nrow(T2)] = "Error"
    class(T2) = "anova"
    Result[[1]] = T2
    Result[[2]] = T4
    names(Result) = c("ANOVA", "EMS")
  }

  if (fQ) attr(Result, "Caution") = "Test assumption: Q's in the corresponding EMS row are 0."
  return(Result)
}

