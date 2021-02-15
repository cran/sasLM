\name{NEWS}
\title{News for Package \pkg{sasLM}}

\section{Version 0.5.0 (2020-02-15)}{
  \itemize{
    \item{\code{lr0} function for simple linear regressions with each independent variable}    
  }
}


\section{Version 0.4.1 (2020-02-15)}{
  \itemize{
    \item{\code{lr} function for linear regression with g2 inverse}    
    \item{\code{tsum} warpper function for \code{tsum0}, \code{tsum1}, \code{tsum2}, \code{tsum3}}    
    \item{\code{bt} function removed to reomove dependency on \code{knitr}}
    \item{\code{NOINT} option in REG, ModelMatrix is removed. Use conventional R formula.}
    \item{Order of coefficients (parameters) changed to be more similar to SAS.}
    \item{\code{T3MS} revived.}
    \item{\code{EMS} has different arguments from T3MS.}
  }
}

\section{Version 0.4.0 (2020-01-29)}{
  \itemize{
    \item{\code{GLM} function now has \code{lsm} option for least square mean}
    \item{\code{LSM} function added for Least Square Means in SAS}
    \item{\code{PDIFF} to do pairwise comparison}
    \item{\code{BY} function to run functions by levels of a variable}
    \item{\code{tsum0}, \code{tsum1}, \code{tsum2}, \code{tsum3} functions for the summary table of a continuous dependent variable}
    \item{\code{bk} function to beautify the output of knitr::kable()}
    \item{\code{bt} function to beautify the output of tsum2, tsum3}
    \item{\code{G2SWEEP} function epxplained}
    \item{\code{regD} function optimized}
    \item{\code{T3MS} renamed as \code{EMS}, because it can handle other types also.}
    \item{\code{CIest} function's arguments and return values changed.}
  }
}

\section{Version 0.3.0 (2020-01-16)}{
  \itemize{
    \item{\code{regD} function added for influence diagnostics}
  }
}

\section{Version 0.2.1 (2020-09-27)}{
  \itemize{
    \item{\code{T3test} a bug fixed for one row table}
    \item{\code{CIest} added for confidence interval estimation}
  }
}

\section{Version 0.2.0 (2020-09-18)}{
  \itemize{
    \item{\code{e3} function revised for RANDOM test functionality}
    \item{\code{satt} added for Satterthwaite approximation}
    \item{\code{T3MS} added for Type 3 Expected Mean Square formula}
    \item{\code{T3test} added for hypothesis test with different error term}
  }
}

\section{Version 0.1.4 (2020-06-02)}{
  \itemize{
    \item{Unnecessary lines in source removed}
    \item{Contrast input for \code{cSS} changed from column vectors to row vectors.}
    \item{\code{REG} can return \code{lfit} result with \code{summarize=FALSE} argument.}
  }
}

\section{Version 0.1.3 (2020-05-13)}{
  \itemize{
    \item{Type III SS function improved}
    \item{Miscellaneous errata corrected}
  }
}

\section{Version 0.1.2 (2020-04-13)}{
  \itemize{
    \item{Code refactored}
  }
}

\section{Version 0.1.1 (2020-04-08)}{
  \itemize{
    \item{An erratum is fixed}
  }
}

\section{Version 0.1.0 (2020-04-08)}{
  \itemize{
    \item{First version}
    \item{TO DO: Function to be added - e4, aov4}
    \item{TO DO: Items to be added - R2, R2ADJ - v0.4.1 \code{lr} do this.}
  }
}
