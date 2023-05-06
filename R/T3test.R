T3test = function(Formula, Data, H="", E="", eps=1e-8)
{
  if (E == "" | length(E) > 1) { stop("One error term should be provided.")
  } else { Error = E }

  if (!attr(terms(Formula, data=Data), "response")) {
    stop("Dependent variable should be provided!")
  }

  x = ModelMatrix(Formula, Data, KeepOrder = FALSE)
  y = model.frame(Formula, Data)[, 1]
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

  ToTest = NULL
  for (i in 1:length(Error)) ToTest = union(ToTest, setdiff(rownames(T0[which(T0[, Error[i]] > 0), ]), Error[i]))
  nTest = length(ToTest)

  Result = list()
  if (nTest > 0) {
    v2 = T1[Error, "Mean Sq"]
    df2 = T1[Error, "Df"]
    for (i in 1:nTest) {
      v1 = T1[ToTest[i], "Mean Sq"]
      df1 = T1[ToTest[i], "Df"]
      w2 = T0[ToTest[i], Error]/T0[Error, Error]
      v4 = satt(vars=c(v2, r1$SSE/r1$DFr), dfs=c(df2, r1$DFr), ws=c(w2, 1 - w2))

      DF = c(df1, v4$Df)
      SS = c(T1[ToTest[i], "Sum Sq"], v4$Variance * v4$Df)
      MS = c(v1, v4$Variance)
      Fval = c(v1/v4$Variance, NA)
      Pval = c(1 - pf(Fval[1], df1, v4$Df), NA)

      T3 = cbind(DF, SS, MS, Fval, Pval)
      colnames(T3) = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
      rownames(T3) = c(ToTest[i], "Error")
      attr(T3, "heading") = paste0("Error: ", round(w2, 4), "*", Error, " + ", round(1 - w2, 4), "*MSE")
      class(T3) = "anova"
      Result[[i]] = T3
    }
    RESIDUALS = c(DF=r1$DFr, SS=r1$SSE, MS=r1$SSE/r1$DFr, Fval=NA, Pval=NA)
    names(RESIDUALS) = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    T2a = T1[setdiff(rownames(T1), ToTest),,drop=FALSE]
    T2b = rbind(T2a, RESIDUALS)
    rownames(T2b) = c(rownames(T2a), "RESIDUALS")
    class(T2b) = "anova"
    Result[[nTest + 1]] = T2b
    names(Result) = c(ToTest, "The Rest Terms")
    if (H != "") {
      Result = Result[H]
      names(Result) = H
    }    
    return(Result)
  } else {
    return(sumANOVA(r1, T1, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x), "factors"))[1]))
  }
}

