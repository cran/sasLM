bk = function(ktab, rpltag=c("n", "N"), dig=10)
{
  if (is.null(rpltag)) {
    for (i in 3:length(ktab)) {
      for (j in dig:1) {
        p0 = paste0(c(".", rep("0", j)), collapse="")
        r0 = paste0(rep(" ", j + 1), collapse="")
        ktab[i] = gsub(p0, r0, ktab[i], fixed=TRUE)
      }
    }
  } else {
    rpltag = union(c("n", "N"), rpltag)
    nRpl = length(rpltag)
    for (i in 1:nRpl) {
      cRpl = rpltag[i]
      ns1 = grep(paste0("[\\|][ ]*", cRpl, "[ ]*[\\|]"), ktab)
      if (length(ns1) > 0) {
        for (j in 1:length(ns1)) {
          for (k in dig:1) {
            p0 = paste0(c(".", rep("0", k)), collapse="")
            r0 = paste0(rep(" ", k + 1), collapse="")
            ktab[ns1[j]] = gsub(p0, r0, ktab[ns1[j]], fixed=TRUE)
          }
        }
      }
    }
  }
  return(ktab)
}
