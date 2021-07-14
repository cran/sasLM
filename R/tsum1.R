tsum1 = function(d, y, u, e=c("Mean", "SD", "N"), ou="", repl=list(c("length"), ("n")))
{
  d = d[!is.na(d[,y]), ]

  if (ou[1] == "") { cNames = as.character(unique(d[, u]))
  } else { cNames = ou }

  nc = length(cNames)
  ne = length(e)

  rNames = e
  for (k in 1:length(repl[[1]])) rNames[rNames == repl[[1]][k]] = repl[[2]][k]

  Res = matrix(nrow=ne, ncol=nc + 1)
  colnames(Res) = c(cNames, "Combined")
  Res = as.data.frame(Res)
  rownames(Res) = rNames

  for (j in 1:nc) {
    tV = d[d[,u] == cNames[j], y]
    for (k in 1:length(e)) Res[k, j] = do.call(e[k], list(tV))
  }

  for (k in 1:length(e)) Res[k, 1 + nc] = do.call(e[k], list(d[,y]))

  attr(Res, "y") = y
  attr(Res, "u") = u
  return(Res)
}

