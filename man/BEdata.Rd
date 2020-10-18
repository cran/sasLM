\name{BEdata}
\alias{BEdata}
\docType{data}
\title{
  An Example Data of Bioequivalence Study
}
\description{
  Contains Cmax data from a real bioequivalence study. 
}
\usage{BEdata}
\format{
  A data frame with 91 observations on the following 6 variables.
  \describe{
    \item{\code{ADM}}{Admission or Hospitalization Group Code: 1, 2, or 3}
    \item{\code{SEQ}}{Group or Sequence character code: 'RT' or 'TR"}
    \item{\code{PRD}}{Period numeric value: 1 or 2}
    \item{\code{TRT}}{Treatment or Drug code: 'R' or 'T'}
    \item{\code{SUBJ}}{Subject ID}
    \item{\code{CMAX}}{Cmax values}
  }
}
\details{
  This contains a real data of 2x2 bioequivalence study, which have three different hospitalization groups. See Bae KS, Kang SH. Bioequivalence data analysis for the case of separate hospitalization. Transl Clin Pharmacol. 2017;25(2):93-100. doi.org/10.12793/tcp.2017.25.2.93
}
\keyword{datasets}
