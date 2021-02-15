tsum3 = function(d, y, l, u, e=c("mean", "sd", "length"), h=NULL, ol1="", ol2="", ou="", rm.dup=TRUE, repl=list(c("length"), c("n")))
{
  if (length(d)*length(y)*length(l)*length(u)*length(e)==0) stop("Check arugments!")
  if (length(l) != 2) stop("length(l) should be 2!")
  if (length(u) != 1) stop("length(u) should be 1!")

  if (is.null(h)) h = list(e, e)

  d = d[!is.na(d[,y]),] # remove NA

  if (ol1[1] == "") { rNames1 = as.character(unique(d[, l[1]]))
  } else { rNames1 = ol1 }
  if (ol2[1] == "") { rNames2 = as.character(unique(d[, l[2]]))
  } else { rNames2 = ol2 }
  if (ou[1] == "") { cNames = as.character(unique(d[, u]))
  } else { cNames = ou }

  nr1 = length(rNames1)
  nr2 = length(rNames2)

  cNames = as.character(unique(d[, u]))
  nc = length(cNames)
  ne = length(e)
  nh1 = length(h[[1]])
  nh2 = length(h[[2]])
  nr = nr1*nr2*ne + nr1*length(h[[1]]) + length(h[[2]])

  Res = matrix(nrow=nr, ncol=nc + 4)
  colnames(Res) = c(l, u, cNames, "Combined")
  Res = as.data.frame(Res)

# Cell Stat
  for (j in 1:nc) {
    for (i1 in 1:nr1) {
      for (i2 in 1:nr2) {
        tV = d[d[,l[1]] == rNames1[i1] & d[,l[2]] == rNames2[i2] & d[,u] == cNames[j], y]

        for (k in 1:ne) {
          cr = (i1 - 1)*(nr2*ne + nh1) + (i2 - 1)*ne + k

          if (rm.dup == TRUE & k > 1) {
            Res[cr, 1] = ""
            Res[cr, 2] = ""
          } else {
            if (i2 == 1) { Res[cr, 1] = rNames1[i1]
            } else { Res[cr, 1] = "" }
            Res[cr, 2] = rNames2[i2]
          }

          Res[cr, 3] = e[k]
          Res[cr, 3 + j] = do.call(e[k], list(tV))
        }
      }

# SubGroup Stat
      tV1 = d[d[,l[1]] == rNames1[i1] & d[,u] == cNames[j], y]
      for (k in 1:nh1) {
        cr = i1*nr2*ne + (i1 - 1)*nh1 + k
        if (rm.dup == TRUE & k > 1) {
          Res[cr, 1] = ""
          Res[cr, 2] = ""
        } else {
          if (rm.dup == TRUE & i2 > 1) { Res[cr, 1] = ""
          } else { Res[cr, 1] = rNames1[i1] }
          Res[cr, 2] = "Combined"
        }
        Res[cr, 3] = h[[1]][k]
        Res[cr, 3 + j] = do.call(h[[1]][k], list(tV1))
      }
    }
  }

# Last Rows Block
  for (j in 1:nc) {
    tV2 = d[d[,u] == cNames[j], y]
    for (k in 1:nh2) {
      cr = nr1*(nr2*ne + nh1) + k
      if (rm.dup == TRUE & k > 1) {
        Res[cr, 1] = ""
        Res[cr, 2] = ""
      } else {
        Res[cr, 1] = "Combined"
        Res[cr, 2] = ""
      }
      Res[cr, 3] = h[[2]][k]
      Res[cr, 3 + j] = do.call(h[[2]][k], list(tV2))
    }
  }

# Last Column except Last Block
  for (i1 in 1:nr1) {
    for (i2 in 1:nr2) {
      tV3 =  d[d[,l[1]] == rNames1[i1] & d[,l[2]] == rNames2[i2], y]
      for (k in 1:ne) {
        cr = (i1 - 1)*(nr2*ne + nh1) + (i2 - 1)*ne + k
        Res[cr, 4 + nc] = do.call(e[k], list(tV3))
      }
    }

    tV4 =  d[d[,l[1]] == rNames1[i1], y]
    for (k in 1:nh1) {
      cr = i1*nr2*ne + (i1 - 1)*nh1 + k
      Res[cr, 4 + nc] = do.call(h[[1]][k], list(tV4))
    }
  }

# Last Block
  for (k in 1:nh2) {
    cr = nr1*nr2*ne + nr1*nh1 + k
    Res[cr, 4 + nc] = do.call(h[[2]][k], list(d[,y]))
  }

# Relpace characters
  for (j in 1:3) {
    for (k in 1:length(repl[[1]])) {
      Res[Res[,j]==repl[[1]][k], j] = repl[[2]][k]
    }
  }

  return(Res)
}

