# IQLM.R
# Installation Qualification (IQ) report generator for the 'sasLM' package.
# Produces a self-contained PDF certifying that sasLM is correctly installed,
# intact, loadable, and operational in the user's own R environment. Uses only base
# R and the package's own pdf helpers (QualReportLM.R), so it needs NO LaTeX/pandoc.

#-- internal: relative-difference comparator ------------------------------------
.iqLMClose = function(observed, reference, tol=1e-3)
{
  observed = as.numeric(observed); reference = as.numeric(reference)
  if (is.na(observed) || is.na(reference)) return(FALSE)
  if (reference == 0) return(abs(observed) <= tol)
  return(abs((observed - reference)/reference) <= tol)
}

#-- internal: parse a DESCRIPTION Depends/Imports string into name/op/version ---
.iqLMParseDeps = function(str)
{
  if (is.null(str) || is.na(str) || !nzchar(str)) return(list())
  parts = trimws(strsplit(str, ",")[[1]])
  out = list()
  for (p in parts) {
    m = regmatches(p, regexec("^([A-Za-z][A-Za-z0-9.]*)[[:space:]]*(\\([[:space:]]*([<>=!]+)[[:space:]]*([0-9.-]+)[[:space:]]*\\))?", p))[[1]]
    if (length(m) >= 2 && nzchar(m[2]))
      out[[length(out) + 1]] = list(name=m[2],
                                    op=if (length(m) >= 4) m[4] else "",
                                    ver=if (length(m) >= 5) m[5] else "")
  }
  return(out)
}

#-- main ------------------------------------------------------------------------
IQLM = function(fileName="sasLM-IQ-Report.pdf", pkgs="sasLM",
                functional=TRUE, performedBy="", paper="auto", sigField=FALSE)
{
  TOL = 1e-3
  now = Sys.time(); si = Sys.info()

  chk = data.frame(Section=character(0), Item=character(0),
                   Result=character(0), Status=character(0), stringsAsFactors=FALSE)
  add = function(section, item, result, status="INFO")
    chk[nrow(chk) + 1, ] <<- list(section, item, as.character(result), status)

  expectedExports = list(
    sasLM = c("GLM", "REG", "aov1", "aov3", "TTEST", "UNIV"))

  ## 1. Test environment ------------------------------------------------------
  S = "1. Test Environment"
  add(S, "R version", R.version.string)
  add(S, "Platform", R.version$platform)
  add(S, "Operating system", paste(si["sysname"], si["release"]))
  add(S, "Machine / node", paste(si["machine"], "/", si["nodename"]))
  add(S, "User (login)", si["login"])
  add(S, "Locale (LC_CTYPE)", Sys.getlocale("LC_CTYPE"))
  libs = .libPaths(); for (i in seq_along(libs)) add(S, paste0("Library path [", i, "]"), libs[i])

  ## 2. Package installation & versions --------------------------------------
  S = "2. Package Installation and Versions"
  allOk = TRUE
  basePkgs = c("base","utils","graphics","grDevices","stats","methods","datasets",
               "tools","grid","splines","stats4","tcltk","compiler","parallel")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly=TRUE)) {
      add(S, paste0(p, ": installed"), "NOT FOUND", "FAIL"); allOk = FALSE; next
    }
    add(S, paste0(p, ": version"),
        tryCatch(as.character(packageVersion(p)), error=function(e) "unknown"), "PASS")
    add(S, paste0(p, ": location"),
        tryCatch(find.package(p), error=function(e) "unknown"), "INFO")
    deps = c(.iqLMParseDeps(packageDescription(p)$Depends),
             .iqLMParseDeps(packageDescription(p)$Imports))
    for (d in deps) {
      if (identical(d$name, "R")) {
        if (!nzchar(d$ver)) next
        ok = do.call(d$op, list(getRversion(), d$ver))
        add(S, paste0(p, " requires R ", d$op, " ", d$ver),
            paste("found", as.character(getRversion())), if (isTRUE(ok)) "PASS" else "FAIL")
        if (!isTRUE(ok)) allOk = FALSE; next
      }
      if (d$name %in% basePkgs) next
      if (!requireNamespace(d$name, quietly=TRUE)) {
        add(S, paste0(p, " requires ", d$name), "dependency NOT installed", "FAIL")
        allOk = FALSE; next
      }
      dv = packageVersion(d$name)
      if (nzchar(d$ver)) {
        ok = do.call(d$op, list(dv, d$ver))
        add(S, paste0(p, " requires ", d$name, " ", d$op, " ", d$ver),
            paste("found", as.character(dv)), if (isTRUE(ok)) "PASS" else "FAIL")
        if (!isTRUE(ok)) allOk = FALSE
      } else add(S, paste0(p, " depends on ", d$name), paste("found", as.character(dv)), "PASS")
    }
  }

  ## 3. File integrity (MD5 checksums) ---------------------------------------
  S = "3. File Integrity (MD5 checksums)"
  for (p in pkgs) {
    if (!requireNamespace(p, quietly=TRUE)) next
    dirp = find.package(p)
    nfile = length(setdiff(list.files(dirp, recursive=TRUE, all.files=TRUE, no.. = TRUE), "MD5"))
    md5 = tryCatch(checkMD5sums(p, dirp), error=function(e) NA)
    if (isTRUE(md5)) {
      add(S, paste0(p, ": checkMD5sums"),
          sprintf("all %d installed files match the MD5 manifest", nfile), "PASS")
    } else if (isTRUE(!md5)) {
      add(S, paste0(p, ": checkMD5sums"),
          "MISMATCH - one or more files differ from the MD5 manifest", "FAIL"); allOk = FALSE
    } else {
      add(S, paste0(p, ": checkMD5sums"),
          paste0("no MD5 manifest in this install - run sasLM::writeMD5LM(\"", p,
                 "\") once after installation (CRAN installs include it)"), "WARN")
    }
  }

  ## 4. Load test and exported functions -------------------------------------
  S = "4. Namespace Load and Exported Functions"
  for (p in pkgs) {
    loaded = tryCatch(requireNamespace(p, quietly=TRUE), error=function(e) FALSE)
    add(S, paste0(p, ": namespace loads"), if (loaded) "OK" else "FAILED",
        if (loaded) "PASS" else "FAIL")
    if (!loaded) { allOk = FALSE; next }
    exps = getNamespaceExports(p); want = expectedExports[[p]]
    if (is.null(want)) {
      add(S, paste0(p, ": exported objects"), paste(length(exps), "objects exported"),
          if (length(exps) > 0) "PASS" else "FAIL")
    } else {
      missing = setdiff(want, exps)
      if (length(missing) == 0) add(S, paste0(p, ": core exports present"),
                                    paste(want, collapse=", "), "PASS")
      else { add(S, paste0(p, ": core exports present"),
                 paste("MISSING:", paste(missing, collapse=", ")), "FAIL"); allOk = FALSE }
    }
  }

  ## 5. Functional verification ----------------------------------------------
  if (isTRUE(functional) && requireNamespace("sasLM", quietly=TRUE)) {
    S = "5. Functional Verification"
    # Goodnight (1976) p7, balanced 2x2 factorial. Reference = SAS PROC GLM,
    # as documented in Validation-Report-GLM (sec 4). External to this computation.
    d = data.frame(A=factor(c(1,1,1,1,2,2,2,2)),
                   B=factor(c(1,1,2,2,1,1,2,2)),
                   y=c(4.91,2.04,5.87,4.42,5.89,7.43,6.84,6.38))
    g = tryCatch(GLM(y ~ A + B + A:B, d), error=function(e) NULL)
    if (is.null(g)) { add(S, "GLM smoke test (Goodnight p7)", "GLM() raised an error", "FAIL"); allOk = FALSE }
    else {
      # extract defensively: a renamed/absent row or column yields NA (-> FAIL via
      # .iqLMClose), so a future regression degrades to NOT QUALIFIED, never a crash.
      pick = function(tab, r, c)
        tryCatch(suppressWarnings(as.numeric(g[[tab]][r, c])), error=function(e) NA_real_)
      checks = list(
        list("ANOVA MODEL Sum Sq",        pick("ANOVA", "MODEL", "Sum Sq"),     13.60265),
        list("ANOVA RESIDUALS Sum Sq",    pick("ANOVA", "RESIDUALS", "Sum Sq"),  6.46130),
        list("Type I  A Sum Sq",          pick("Type I", "A", "Sum Sq"),        10.81125),
        list("Type III A:B Sum Sq",       pick("Type III", "A:B", "Sum Sq"),     1.47920),
        list("Fitness R-square",          pick("Fitness", 1, "R-square"),        0.6779647))
      for (ck in checks) {
        obs = suppressWarnings(as.numeric(ck[[2]])); ok = .iqLMClose(obs, ck[[3]], TOL)
        add(S, ck[[1]], sprintf("obs=%.6g  ref=%.6g", obs, ck[[3]]),
            if (ok) "PASS" else "FAIL")
        if (!ok) allOk = FALSE
      }
    }
    c1 = .iqLMClose(1.234567, 1.234567, TOL)
    add(S, "Control 1: identical values compare equal",
        if (c1) "comparator returned EQUAL" else "comparator failed",
        if (isTRUE(c1)) "PASS" else "FAIL"); if (!isTRUE(c1)) allOk = FALSE
    c2 = .iqLMClose(1.00, 1.10, TOL)
    add(S, "Control 2: 10% different values flagged",
        if (!c2) "comparator returned DIFFERENT" else "comparator failed",
        if (isFALSE(c2)) "PASS" else "FAIL"); if (!isFALSE(c2)) allOk = FALSE
  }

  ## ---- assemble report text (page 1 = signatures, then detail pages) --------
  hr  = paste(rep("=", 72), collapse=""); hr2 = paste(rep("-", 72), collapse="")
  pkgVer = function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "not installed")
  nFail = sum(chk$Status == "FAIL"); nWarn = sum(chk$Status == "WARN"); nPass = sum(chk$Status == "PASS")
  qualified = isTRUE(allOk) && nFail == 0
  status = if (qualified) "QUALIFIED" else "NOT QUALIFIED"
  who = if (nzchar(performedBy)) performedBy else as.character(si["login"])
  stamp = format(now, "%Y-%m-%d %H:%M:%S")

  sig = character(0); sg = function(...) sig <<- c(sig, ...)
  sg(hr); sg("        INSTALLATION QUALIFICATION (IQ) REPORT"); sg(hr); sg("")
  sg("Qualified package(s):")
  for (p in pkgs) sg(sprintf("    %-14s %s", p, pkgVer(p)))
  sg("")
  sg(sprintf("Report generated : %s %s", stamp, Sys.timezone()))
  sg(sprintf("Generated by     : sasLM::IQLM(), sasLM %s", pkgVer("sasLM")))
  sg(sprintf("Performed by     : %s", who))
  sg("")
  sg(sprintf("Overall result   : %s", status))
  sg(sprintf("                   %d PASS, %d FAIL, %d WARN", nPass, nFail, nWarn))
  sg("")
  if (qualified)
    sg(strwrap(paste("Based on the checks documented on the following pages,",
                     "sasLM is correctly installed, intact, loadable, and",
                     "operational in this R environment as of the date above."), width=72))
  else
    sg(strwrap(paste("One or more checks failed (see the following pages). The",
                     "installation is NOT qualified until the failures are resolved",
                     "and this report is re-generated."), width=72))
  if (nWarn > 0) {
    sg("")
    sg(strwrap(paste("Note: [WARN] items did not fail qualification but should be",
                     "reviewed (e.g. a missing MD5 manifest means file integrity",
                     "could not be cryptographically verified for that package)."), width=72))
  }
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
  sg(""); sg(""); sg("")
  sg("Detailed qualification evidence follows on the subsequent pages.")

  body = character(0); ad = function(...) body <<- c(body, ...)
  wrapVal = function(prefix, value, contIndent) {
    width = 72; avail = width - nchar(prefix); if (avail < 8) avail = 8
    if (nchar(value) <= avail) return(paste0(prefix, value))
    chunks = character(0)
    while (nchar(value) > 0) {
      take = min(avail, nchar(value)); chunks = c(chunks, substr(value, 1, take))
      value = substr(value, take + 1, nchar(value)); avail = width - contIndent
    }
    pad = paste(rep(" ", contIndent), collapse="")
    c(paste0(prefix, chunks[1]), if (length(chunks) > 1) paste0(pad, chunks[-1]))
  }
  for (S in unique(chk$Section)) {
    ad(""); ad(S); ad(hr2)
    sub = chk[chk$Section == S, ]
    for (i in seq_len(nrow(sub))) {
      prefix = sprintf(" [%-4s] %-36s : ", sub$Status[i], substr(sub$Item[i], 1, 36))
      for (ln in wrapVal(prefix, sub$Result[i], nchar(prefix))) ad(ln)
    }
  }
  ad("\f"); ad(hr); ad("        OVERALL RESULT"); ad(hr)
  ad(sprintf("    STATUS: %s   (%d PASS, %d FAIL, %d WARN)", status, nPass, nFail, nWarn))
  ad("\f"); ad("Appendix A. sessionInfo()"); ad(hr2)
  ow = options(width=78)
  siLines = tryCatch(capture.output(print(sessionInfo())),
                     error=function(e) "sessionInfo() unavailable")
  options(ow)
  for (ln in siLines) ad(ln)
  ad("\f"); ad("Appendix B. Installed File Checksums (md5)"); ad(hr2)
  ad("md5sum of every installed file, computed in this environment.")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly=TRUE)) next
    dirp = find.package(p)
    files = setdiff(list.files(dirp, recursive=TRUE, all.files=TRUE, no.. = TRUE), "MD5")
    sums = tryCatch(md5sum(file.path(dirp, files)), error=function(e) rep(NA, length(files)))
    ad(""); ad(sprintf("%s %s  (%d files)", p, pkgVer(p), length(files)))
    pth = paste0("at ", dirp)
    while (nchar(pth) > 0) { ad(substr(pth, 1, 72)); pth = substr(pth, 73, nchar(pth)) }
    ad(paste(rep("-", 72), collapse=""))
    for (i in order(files)) {
      s = if (is.na(sums[i])) paste(rep("?", 32), collapse="") else sums[i]
      f = files[i]; if (nchar(f) > 38) f = paste0("...", substr(f, nchar(f) - 34, nchar(f)))
      ad(sprintf("%s  %s", s, f))
    }
  }

  rinfo = .qRenderPDF(fileName, sigLines=sig, bodyLines=body,
                      title="Installation Qualification Report", verdict=status,
                      footer=sprintf("sasLM::IQLM  %s", stamp), paper=paper)
  if (isTRUE(sigField))
    try(addSigFieldLM(fileName, page=1L, fieldNames=c("Performed_by", "Reviewed_by"),
                      rects=list(.qSigRect(rinfo, perfRow), .qSigRect(rinfo, revRow))),
        silent=TRUE)
  apath = normalizePath(fileName, winslash="/", mustWork=FALSE)
  cat(sprintf("IQ report (STATUS: %s, paper: %s) saved at:\n  %s\n",
              status, rinfo$paper, apath))
  if (isTRUE(sigField))
    cat("  Open this PDF in Adobe Acrobat Reader and click a signature field to sign.\n")
  invisible(list(fileName=apath, qualified=qualified, paper=rinfo$paper,
                 checks=chk, nPass=nPass, nFail=nFail, nWarn=nWarn))
}
