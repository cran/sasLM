\name{IQLM}
\alias{IQLM}
\title{Installation Qualification (IQ) report to a pdf file}
\description{
Generate a self-contained PDF report documenting whether the \pkg{sasLM} package is
correctly installed, intact, loadable, and operational in the user's own R
environment. It is intended as Installation Qualification (IQ) evidence, in the
spirit of commercial statistical-software validation. The report uses only base R and
the package's own pdf helpers, so it requires no LaTeX, pandoc, or other external
tools. The signature (approval) page is placed first and the report uses 1 inch
margins on every side, with Letter paper in a United States locale and A4 elsewhere.
}
\usage{
IQLM(fileName = "sasLM-IQ-Report.pdf", pkgs = "sasLM",
     functional = TRUE, performedBy = "", paper = "auto", sigField = FALSE)
}
\arguments{
  \item{fileName}{file name to save the PDF report.}
  \item{pkgs}{character vector of package names to qualify; defaults to
    \code{"sasLM"}. Declared dependencies are checked automatically.}
  \item{functional}{if \code{TRUE}, run a small operational check: a known linear
    model (Goodnight 1976 p7, a balanced 2x2 factorial) is fitted with
    \code{\link{GLM}} and selected results (model and residual sums of squares, a
    Type I and a Type III sum of squares, and R-square) are compared with the SAS
    PROC GLM reference within a relative tolerance of 1e-3, plus two control tests
    that verify the comparator itself.}
  \item{performedBy}{name of the person performing the qualification, printed on the
    signature page. Defaults to the login name.}
  \item{paper}{paper size: \code{"letter"}, \code{"a4"}, or \code{"auto"} (the
    default; Letter in a United States locale, A4 elsewhere).}
  \item{sigField}{if \code{TRUE}, add Adobe Acrobat Reader signature fields to the
    finished report (via \code{\link{addSigFieldLM}}) so it can be signed with one
    click in Acrobat Reader.}
}
\details{
The report contains: the test environment (R version, platform, OS, locale, library
paths); the installed package version, location and declared-dependency satisfaction;
file integrity via \code{\link[tools]{checkMD5sums}} (PASS / FAIL / WARN when no
manifest is present - see \code{\link{writeMD5LM}}); namespace load and core exports;
an optional functional verification; an overall QUALIFIED / NOT QUALIFIED verdict; a
\code{\link[utils]{sessionInfo}} appendix; and a per-file md5 checksum appendix.
}
\value{
Invisibly, a list with \code{fileName}, \code{qualified}, \code{paper}, \code{checks}
(a data frame of every check), and the counts \code{nPass}, \code{nFail},
\code{nWarn}. A PDF file is written to \code{fileName}.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{OQLM}}, \code{\link{writeMD5LM}}, \code{\link{GLM}}
}
\examples{
#IQLM()
#res <- IQLM(performedBy = "Jane Doe"); res$qualified
}
\keyword{misc}
