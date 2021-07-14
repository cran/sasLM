tsum2 = function(d, y, l, u, e=c("Mean", "SD", "N"), h=NULL, ol="", ou="", rm.dup=TRUE, repl=list(c("length"), c("n")))
{
  if (length(d)*length(y)*length(l)*length(u)*length(e)==0) stop("Check arugments!")
  if (length(l) != 1) stop("length(l) should be 1!")
  if (length(u) != 1) stop("length(u) should be 1!")

  if (is.null(h)) h = e

  d = d[!is.na(d[,y]),] # remove NA

  if (ol[1] == "") { rNames = as.character(unique(d[, l]))
  } else { rNames = ol }
  if (ou[1] == "") { cNames = as.character(unique(d[, u]))
  } else { cNames = ou }

  nr = length(rNames)
  nc = length(cNames)
  ne = length(e)

  Res = matrix(nrow=nr*ne + length(h), ncol=nc + 3)
  colnames(Res) = c(l, u, cNames, "Combined")
  Res = as.data.frame(Res)

# Cell Stat
  for (j in 1:nc) {
    for (i in 1:nr) {
      tV = d[d[,l] == rNames[i] & d[,u] == cNames[j], y]
      for (k in 1:ne) {
        cr = (i - 1)*ne + k
        if (rm.dup == TRUE & k > 1) {
          Res[cr, 1] = ""
        } else {
          Res[cr, 1] = rNames[i]
        }
        Res[cr, 2] = e[k]
        Res[cr, 2 + j] = do.call(e[k], list(tV))
      }
    }
# Group Stat
    tV1 = d[d[,u] == cNames[j], y]
    for (k in 1:length(h)) {
      cr = nr*ne + k
      if (rm.dup == TRUE & k > 1) {
        Res[cr, 1] = ""
      } else {
        Res[cr, 1] = "Combined"
      }
      Res[cr, 2] = h[k]
      Res[cr, 2 + j] = do.call(h[k], list(tV1))
    }
  }

# Last Column except Last Block
  for (i in 1:nr) {
    tV2 = d[d[,l] == rNames[i], y]
    for (k in 1:ne) {
      cr = i*ne - ne + k
      Res[cr, nc+3] = do.call(e[k], list(tV2))
    }
  }

# Last Block
  for (k in 1:length(h)) {
    Res[nr*ne + k, 3 + nc] = do.call(h[k], list(d[, y]))
  }

# Relpace characters
  for (j in 1:2) {
    for (k in 1:length(repl[[1]])) {
      Res[Res[,j]==repl[[1]][k], j] = repl[[2]][k]
    }
  }

  return(Res)
}
