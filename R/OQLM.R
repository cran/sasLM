# OQLM.R
# Operational Qualification (OQ) report generator for the 'sasLM' linear-model
# engine. Runs a set of GLM test scenarios on the user's machine and compares
# every computed value (ANOVA table, fit statistics, and Type I/II/III sums of
# squares) to a pre-specified reference value within a relative tolerance, then
# writes a self-contained PDF that lists every value next to its reference. Uses
# only base R + the package's pdf helpers (QualReportLM.R): no LaTeX.
#
# Reference values are the SAS PROC GLM / textbook results documented in the
# package's validation report (Validation-Report-GLM), shipped in inst/OQ as
# per-scenario <id>.ref.csv files with the scenario manifest OQcfg.csv. The first
# release covers GLM(); other procedures can be added by extending inst/OQ.

#-- internal: canonicalize a name for robust matching (spacing/case tolerant).
# ':' is kept so an interaction term (e.g. "A:B") cannot collide with a main
# effect literally named "AB", nor "Sire:Ration" with a factor "SireRation".
.oqLMCanon = function(x) tolower(gsub("[^A-Za-z0-9:]", "", x))

#-- internal: symmetric relative difference, with absolute fallback near 0 -------
.oqLMDiff = function(ref, obs, tol, absTol=1e-6)
{
  ref = suppressWarnings(as.numeric(ref)); obs = suppressWarnings(as.numeric(obs))
  if (is.na(ref) && is.na(obs)) return(list(rd=0, pass=TRUE))
  if (is.na(ref) ||  is.na(obs)) return(list(rd=NA_real_, pass=FALSE))
  denom = (abs(ref) + abs(obs)) / 2
  if (denom <= absTol) return(list(rd=0, pass=(abs(ref - obs) <= absTol)))  # near-zero target: absolute check
  rd = abs(ref - obs) / denom
  list(rd=rd, pass=(rd <= tol))
}

#-- internal: flatten a GLM result to long (table, term, stat, value) ----------
.oqLMFlatten = function(g)
{
  tabmap = c(ANOVA="ANOVA", Fitness="Fitness",
             "Type I"="TypeI", "Type II"="TypeII", "Type III"="TypeIII")
  tab = character(0); term = character(0); stat = character(0); val = numeric(0)
  for (nm in names(g)) {
    key = tabmap[nm]; if (is.na(key)) next
    M = g[[nm]]
    if (key == "Fitness") {
      v = suppressWarnings(as.numeric(M[1, ])); cn = colnames(M)
      tab = c(tab, rep(key, length(v))); term = c(term, rep("", length(v)))
      stat = c(stat, cn); val = c(val, v)
    } else {
      Mm = as.matrix(M); rn = rownames(Mm); cn = colnames(Mm)
      for (r in rn) for (cc in cn) {
        tab = c(tab, key); term = c(term, r); stat = c(stat, cc)
        val = c(val, suppressWarnings(as.numeric(Mm[r, cc])))
      }
    }
  }
  data.frame(table=tab, term=term, stat=stat, value=val,
             key=paste(.oqLMCanon(tab), .oqLMCanon(term), .oqLMCanon(stat), sep="|"),
             stringsAsFactors=FALSE)
}

#-- internal: run a single scenario -> full comparison table --------------------
.oqLMRunCase = function(cfgRow, refDir, tol)
{
  id = cfgRow$id
  out = list(id=id, desc=cfgRow$desc, source=cfgRow$source, status="FAIL",
             nComp=0, nFail=0, maxRD=NA_real_, comp=NULL, settings="", message="")
  out$settings = sprintf("data=%s, formula=%s, factors=%s",
                         cfgRow$datafile, cfgRow$formula, cfgRow$factors)

  d = tryCatch(read.csv(file.path(refDir, cfgRow$datafile), check.names=TRUE),
               error=function(e) NULL)
  if (is.null(d)) { out$message = paste("data not found:", cfgRow$datafile); return(out) }
  facs = trimws(strsplit(cfgRow$factors, ";")[[1]]); facs = facs[nzchar(facs)]
  d = tryCatch(af(d, facs), error=function(e) d)

  g = tryCatch(GLM(stats::as.formula(cfgRow$formula), d), error=function(e) e)
  if (inherits(g, "error")) { out$message = paste("GLM error:", conditionMessage(g)); return(out) }

  W = tryCatch(read.csv(file.path(refDir, cfgRow$reffile), as.is=TRUE, check.names=FALSE),
               error=function(e) NULL)
  if (is.null(W)) { out$message = paste("reference not found:", cfgRow$reffile); return(out) }

  flat = .oqLMFlatten(g)
  if (anyDuplicated(flat$key)) {     # defensive: ambiguous labels would mis-match silently
    out$message = "ambiguous result labels (duplicate canonical keys); cannot match reliably"
    return(out)
  }
  Wkey = paste(.oqLMCanon(W$table), .oqLMCanon(W$term), .oqLMCanon(W$stat), sep="|")
  nR = nrow(W)
  obs = rep(NA_real_, nR); rd = rep(NA_real_, nR); ps = logical(nR)
  for (i in seq_len(nR)) {
    j = match(Wkey[i], flat$key)
    o = if (is.na(j)) NA_real_ else flat$value[j]
    obs[i] = o
    st = .oqLMDiff(W$value[i], o, tol)
    rd[i] = st$rd; ps[i] = isTRUE(st$pass)
  }
  comp = data.frame(table=W$table, term=W$term, stat=W$stat,
                    ref=suppressWarnings(as.numeric(W$value)), obs=obs,
                    reldiff=rd, pass=ps, stringsAsFactors=FALSE)
  out$comp = comp
  out$nComp = nR; out$nFail = sum(!ps)
  out$maxRD = if (all(is.na(rd))) NA_real_ else max(rd, na.rm=TRUE)
  out$status = if (nR > 0 && out$nFail == 0) "PASS" else "FAIL"
  out
}

#-- main ------------------------------------------------------------------------
OQLM = function(fileName="sasLM-OQ-Report.pdf", cfg=NULL, tol=1e-3,
                refDir=system.file("OQ", package="sasLM"), performedBy="",
                paper="auto", sigField=FALSE)
{
  now = Sys.time(); si = Sys.info()
  if (!nzchar(refDir) || !dir.exists(refDir))
    stop("OQ reference directory not found. Install 'sasLM' or pass refDir=.")
  if (is.null(cfg)) {
    cfg = tryCatch(read.csv(file.path(refDir, "OQcfg.csv"), as.is=TRUE), error=function(e) NULL)
    if (is.null(cfg)) stop("OQcfg.csv not found in refDir.")
  }
  needCols = c("id","desc","datafile","formula","factors","reffile","source")
  if (!all(needCols %in% colnames(cfg)))
    stop("OQcfg must have columns: ", paste(needCols, collapse=", "))

  results = lapply(seq_len(nrow(cfg)), function(i) .oqLMRunCase(cfg[i, ], refDir, tol))

  pkgVer = function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "not installed")
  hr  = paste(rep("=", 72), collapse=""); hr2 = paste(rep("-", 72), collapse="")
  nFailCases = sum(vapply(results, function(r) r$status != "PASS", logical(1)))
  qualified  = (nFailCases == 0) && length(results) > 0
  status = if (qualified) "QUALIFIED" else "NOT QUALIFIED"
  who = if (nzchar(performedBy)) performedBy else as.character(si["login"])
  stamp = format(now, "%Y-%m-%d %H:%M:%S")
  nCompAll = sum(vapply(results, function(r) r$nComp, integer(1)))
  nFailAll = sum(vapply(results, function(r) r$nFail, integer(1)))

  sig = character(0); sg = function(...) sig <<- c(sig, ...)
  sg(hr); sg("        OPERATIONAL QUALIFICATION (OQ) REPORT"); sg(hr); sg("")
  sg("Engine under test:")
  sg(sprintf("    %-14s %s", "sasLM", pkgVer("sasLM")))
  sg("")
  sg(sprintf("Report generated : %s %s", stamp, Sys.timezone()))
  sg(sprintf("Performed by     : %s", who))
  sg(sprintf("Acceptance       : symmetric relative difference <= %g", tol))
  sg("Reference         : SAS PROC GLM / textbook values (Validation-Report-GLM)")
  sg("")
  sg(sprintf("Overall result   : %s", status))
  sg(sprintf("                   %d of %d scenarios passed (%d of %d values)",
             length(results) - nFailCases, length(results), nCompAll - nFailAll, nCompAll))
  sg("")
  if (qualified)
    sg(strwrap(paste("Every computed value in all operational scenarios reproduced",
                     "the pre-specified reference value within the acceptance",
                     "tolerance on this system. Per-scenario details follow."), width=72))
  else
    sg(strwrap(paste("One or more values did not meet the acceptance criterion.",
                     "See the per-scenario details on the following pages."), width=72))
  sg(""); sg(""); sg("Approval"); sg(hr2)
  sg("Sign this report in Adobe Acrobat Reader (no print/scan): open it,")
  sg("choose \"Use a certificate\" > \"Digitally sign\", sign with your Digital")
  sg("ID, then click the signature box under each name below. (Generate with")
  sg("sigField=TRUE, or run addSigFieldLM(), to add the boxes.) Scriptable")
  sg("alternative: signPDFLM()/verifyPDFLM().")
  sg("")
  perfRow = length(sig) + 1L
  sg("    Performed by : ____________________________   Date: __________")
  sg(""); sg(""); sg("")
  revRow = length(sig) + 1L
  sg("    Reviewed by  : ____________________________   Date: __________")
  sg(""); sg("")

  body = character(0); ad = function(...) body <<- c(body, ...)
  ad("Scenario Summary"); ad(hr2)
  ad(sprintf(" %-4s %-40s %6s %5s %11s", "Stat", "Scenario", "Values", "Fail", "max rel.df")); ad(hr2)
  for (r in results) {
    mrd = if (is.na(r$maxRD)) "  NA" else sprintf("%.2e", r$maxRD)
    ad(sprintf(" [%-2s] %-40s %6d %5d %11s", if (r$status == "PASS") "OK" else "X",
               substr(r$id, 1, 40), r$nComp, r$nFail, mrd))
  }
  ad("")
  ad("Scenario Details - every computed value vs its reference"); ad(hr2)
  fmtNum = function(x) if (is.na(x)) "NA" else formatC(x, format="g", digits=8)
  fmtE   = function(x) if (is.na(x)) "NA" else sprintf("%.2e", x)
  tabLabel = c(ANOVA="ANOVA", Fitness="Fitness", TypeI="Type I", TypeII="Type II", TypeIII="Type III")
  hdr = sprintf("   %-22s %12s %12s %9s %9s %2s",
                "TERM / STATISTIC", "REFERENCE", "OBSERVED", "ABS.DIFF", "REL.DIFF", "ST")
  for (r in results) {
    ad("\f"); ad(hr); ad(sprintf("Scenario: %s   [%s]", r$id, r$status))
    for (ln in strwrap(r$desc, width=72, exdent=3)) ad(ln)
    for (ln in strwrap(paste("Source:", r$source), width=72, exdent=3)) ad(ln)
    for (ln in strwrap(r$settings, width=72, exdent=3)) ad(ln)
    if (nzchar(r$message)) { ad(sprintf("   note: %s", r$message)); next }
    ad(sprintf("   %d values compared; %d failed; max rel.diff %s",
               r$nComp, r$nFail, fmtE(r$maxRD)))
    cp = r$comp
    for (tb in c("ANOVA","Fitness","TypeI","TypeII","TypeIII")) {
      idx = which(cp$table == tb)
      if (length(idx) == 0) next
      ad(""); ad(sprintf(" %s", tabLabel[tb])); ad(hdr)
      ad(paste0("   ", paste(rep("-", 69), collapse="")))
      for (i in idx) {
        lab = trimws(paste(cp$term[i], cp$stat[i]))
        ab = abs(cp$ref[i] - cp$obs[i])
        ad(sprintf("   %-22s %12s %12s %9s %9s %2s",
                   substr(lab, 1, 22), fmtNum(cp$ref[i]), fmtNum(cp$obs[i]),
                   fmtE(ab), fmtE(cp$reldiff[i]), if (cp$pass[i]) "OK" else "X"))
      }
    }
  }
  ad("\f"); ad(hr); ad("        OVERALL RESULT"); ad(hr)
  ad(sprintf("    STATUS: %s", status))
  ad(sprintf("    %d of %d scenarios passed; %d of %d values within tolerance",
             length(results) - nFailCases, length(results), nCompAll - nFailAll, nCompAll))
  ad("\f"); ad("Appendix A. sessionInfo()"); ad(hr2)
  ow = options(width=78)
  siLines = tryCatch(capture.output(print(sessionInfo())), error=function(e) "unavailable")
  options(ow)
  for (ln in siLines) ad(ln)

  rinfo = .qRenderPDF(fileName, sigLines=sig, bodyLines=body,
                      title="Operational Qualification Report", verdict=status,
                      footer=sprintf("sasLM::OQLM  %s", stamp), paper=paper)
  if (isTRUE(sigField))
    try(addSigFieldLM(fileName, page=1L, fieldNames=c("Performed_by", "Reviewed_by"),
                      rects=list(.qSigRect(rinfo, perfRow), .qSigRect(rinfo, revRow))),
        silent=TRUE)
  apath = normalizePath(fileName, winslash="/", mustWork=FALSE)
  cat(sprintf("OQ report (STATUS: %s, paper: %s, %d/%d scenarios passed) saved at:\n  %s\n",
              status, rinfo$paper, length(results) - nFailCases, length(results), apath))
  if (isTRUE(sigField))
    cat("  Open this PDF in Adobe Acrobat Reader and click a signature field to sign.\n")
  invisible(list(fileName=apath, qualified=qualified, tol=tol, results=results,
                 paper=rinfo$paper, nScenarios=length(results), nFailScenarios=nFailCases,
                 nValues=nCompAll, nFailValues=nFailAll))
}
