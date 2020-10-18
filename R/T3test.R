T3test = function(Formula, Data, Error="", eps=1e-8)
{
  if (Error == "") stop("Error term is needed!")
  y = model.frame(Formula, Data)[, 1]
  x = ModelMatrix(Formula, Data, KeepOrder = FALSE)
  r1 = lfit(x, y, eps=eps)

  L0 = e3(Formula, Data, eps=eps)
  T1 = SS(x, r1, L0)

  if (!(Error %in% rownames(T1))) {
    x2 = attr(x$terms, "term.labels")
    xi = grep(":", x2)
    msg1 = "Error term is not found!"
    msg2 = "  Choose an error term among the following: "
    msg3 = paste(x2[xi], collapse=", ")
    stop(paste(msg1, msg2, msg3, sep="\n")) 
  }

  T0 = T3MS(Formula, Data, L0)

  ToTest = setdiff(rownames(T0[which(T0[,Error] > 0),]), Error)
  nTest = length(ToTest)

  Result = list()
  if (nTest > 0) {
    v2 = T1[Error, "Mean Sq"]
    df2 = T1[Error, "Df"]
    for (i in 1:nTest) {
      v1 = T1[ToTest[i], "Mean Sq"]
      df1 = T1[ToTest[i], "Df"]
      w2 = T0[ToTest[i], Error]/T0[Error, Error]
      v4 = satt(c(w2, 1 - w2), c(v2, r1$SSE/r1$DFr), c(df2, r1$DFr))
    
      DF = c(df1, v4$Df)
      SS = c(T1[ToTest[i], "Sum Sq"], v4$Variance * v4$Df)
      MS = c(v1, v4$Variance)
      Fval = c(v1/v4$Variance, NA)
      Pval = c(1 - pf(Fval[1], df1, v4$Df), NA)
    
      T3 = cbind(DF, SS, MS, Fval, Pval)
      colnames(T3) = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
      rownames(T3) = c(ToTest[i], "Error")
      attr(T3, "heading") = paste0("Error: ", format(w2, digits=4), "*", Error, " + ", format(1 - w2, digits=4), "*MSE")       
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
    return(Result)
  } else {
    return(sumANOVA(r1, T1, crossprod(y - mean(y)), nrow(x$X), rownames(attr(terms(x), "factors"))[1]))
  }
}

