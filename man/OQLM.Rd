\name{OQLM}
\alias{OQLM}
\title{Operational Qualification (OQ) report to a pdf file}
\description{
Generate a self-contained PDF report that documents whether the \pkg{sasLM}
linear-model engine reproduces pre-specified reference results on the user's own
machine, within a pre-specified relative tolerance. It is intended as Operational
Qualification (OQ) evidence. For each scenario the report lists every computed value
- the ANOVA table, the fit statistics, and the Type I, II and III sums of squares -
next to its reference value, with the absolute and relative difference and a pass
flag, then an overall verdict. It uses only base R and the package's pdf helpers (no
LaTeX). The signature page is first; 1 inch margins; Letter paper in a United States
locale and A4 elsewhere.
}
\usage{
OQLM(fileName = "sasLM-OQ-Report.pdf", cfg = NULL, tol = 0.001,
     refDir = system.file("OQ", package = "sasLM"), performedBy = "",
     paper = "auto", sigField = FALSE)
}
\arguments{
  \item{fileName}{file name to save the PDF report.}
  \item{cfg}{a data frame of scenario definitions. If \code{NULL}, the manifest
    \file{OQcfg.csv} in \code{refDir} is used. It must have the columns \code{id},
    \code{desc}, \code{datafile} (an input csv in \code{refDir}), \code{formula},
    \code{factors} (semicolon-separated column names made factors with
    \code{\link{af}}), \code{reffile} (the reference csv in \code{refDir}), and
    \code{source} (a citation). Each reference csv has columns \code{table}
    (\code{ANOVA}, \code{Fitness}, \code{TypeI}, \code{TypeII} or \code{TypeIII}),
    \code{term}, \code{stat}, and \code{value}.}
  \item{tol}{acceptance tolerance: a value passes when the symmetric relative
    difference \code{abs(ref - obs) / ((abs(ref) + abs(obs))/2)} is at most
    \code{tol}; an absolute fallback is used when both values are essentially zero.}
  \item{refDir}{directory holding \file{OQcfg.csv}, the input data, and the reference
    csv files. Defaults to the \code{OQ} folder installed with \pkg{sasLM}.}
  \item{performedBy}{name of the person performing the qualification, printed on the
    signature page. Defaults to the login name.}
  \item{paper}{paper size: \code{"letter"}, \code{"a4"}, or \code{"auto"} (the
    default; Letter in a United States locale, A4 elsewhere).}
  \item{sigField}{if \code{TRUE}, add Adobe Acrobat Reader signature fields to the
    finished report (via \code{\link{addSigFieldLM}}) so it can be signed with one
    click in Acrobat Reader.}
}
\details{
For each scenario the function reads the (frozen) input data, makes the named columns
factors with \code{\link{af}}, fits the model with \code{\link{GLM}}, flattens the
result, and compares every reference value to the corresponding computed value. A
scenario passes only when all of its values pass.

The first release covers \code{\link{GLM}} with scenarios from Harvey (ARS20-8) and
Goodnight (1976). The bundled reference values are the SAS PROC GLM / textbook
results documented in the package validation report (\emph{Validation of 'sasLM'
Package}); they establish concordance with SAS, the de-facto reference. A new
\code{GLM} scenario can be added with no code change, by adding an input csv, a
reference csv, and one row to \file{OQcfg.csv} in \code{refDir}. Adding a different
procedure (such as \code{\link{REG}}, \code{\link{aov3}}, \code{\link{TTEST}} or
\code{\link{UNIV}}) additionally requires extending \code{OQLM} itself (procedure
dispatch and result flattening).
}
\value{
Invisibly, a list with \code{fileName}, \code{qualified}, \code{tol},
\code{nScenarios}, \code{nFailScenarios}, \code{nValues}, \code{nFailValues},
\code{paper}, and \code{results} (a per-scenario list whose \code{comp} element is the
full table of every comparison). A PDF is written to \code{fileName}.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{IQLM}}, \code{\link{writeMD5LM}}, \code{\link{GLM}}
}
\examples{
#OQLM()
#res <- OQLM(); res$qualified; res$results[[1]]$comp
}
\keyword{misc}
