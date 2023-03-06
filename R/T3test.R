T3test = function(Formula, Data, Error="", eps=1e-8)
{
  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  x = ModelMatrix(Formula, Data, KeepOrder = FALSE)
  y = model.frame(Formula, Data)[,1]
  if (!is.numeric(y)) stop("Dependent variable should be numeric!")

  r1 = lfit(x, y, eps=eps)
  L0 = e3(x, eps=eps)
  T1 = SS(x, r1, L0)
  T2 = rbind(T1, c(Df=r1$DFr, r1$SSE, r1$SSE/r1$DFr, NA, NA))
  rownames(T2) = c(rownames(T1), "MSE")

  RowNamesT1 = rownames(T1)
  strError = strsplit(Error, ":")
  strTerm = strsplit(RowNamesT1, ":") 
  for (i in 1:nrow(T1)) {
    for (j in 1:length(Error)) {
      if (setequal(strTerm[[i]], strError[[j]])) Error[j] = RowNamesT1[i] 
    }
  }

  if (!all(Error %in% rownames(T1))) {
    x2 = attr(x$terms, "term.labels")
    xi = grep(":", x2)
    msg1 = "Some error terms are not found!"
    msg2 = "  Choose an error term among the following: "
    msg3 = paste(x2[xi], collapse=", ")
    stop(paste(msg1, msg2, msg3, sep="\n"))
  }

  T0 = T3MS(Formula, Data, L0)
  W0 = T0
  diag(W0) = 0
  for (i in 1:nrow(T0)) W0[, i] = W0[, i]/T0[i, i]

## Check nested
  ColNames = strsplit(colnames(W0), ":")
  RowNames = strsplit(rownames(W0), ":")
  for (i in 1:nrow(W0)) {
    for (j in 1:(ncol(W0) - 1)) {
      if (W0[i, j] > 0) {
        for (k in 1:ncol(W0)) {
          if (k == j) next
          if (all(ColNames[[j]] %in% ColNames[[k]])) W0[i, k] = W0[i, k] - 1
        }
      }
    }
  }
##

  W0 = cbind(W0, MSE = 1 - round(apply(W0, 1, sum), 8))

  ToTest = NULL
  for (i in 1:length(Error)) ToTest = union(ToTest, setdiff(rownames(T0[which(T0[, Error[i]] > 0), ]), Error[i]))
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
    return(Result)
  } else {
    return(sumANOVA(r1, T1, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x), "factors"))[1]))
  }
}

