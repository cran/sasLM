\name{NEWS}
\title{News for Package \pkg{sasLM}}

\section{Version 0.9.10 (2023-07-10)}{
  \itemize{
    \item{Typo correction}
  }
}

\section{Version 0.9.9 (2023-06-05)}{
  \itemize{
    \item{TTEST, mtest, tmtest, ztest, vtest functions are added.}
    \item{Every time just before calling mvtnorm::pmvt, mvtnorm::qmvt, set.seed(5) is called.}
  }
}

\section{Version 0.9.8 (2023-05-05)}{
  \itemize{
    \item{RanTest with no fixed effect and only 1 random effect return correct result.}
    \item{T3test supports H and E options}
    \item{LIbin, LIpois, LInorm, and RDLI are removed, and will be included in a separate package.}
  }
}

\section{Version 0.9.7 (2023-04-05)}{
  \itemize{
    \item{Skewness function revised.}
    \item{Bugs of some cases at RDmn1, RDmn, ORmn1, ORmn are fixed.}
    \item{LIbin, LIpois, LInorm, and RDLI are introduced.}
  }
}

\section{Version 0.9.6 (2023-03-06)}{
  \itemize{
    \item{RanTest function added for test with random effect.}
    \item{T3MS will be deprecated soon. Use EMS instead.}
    \item{SLICE function's options are renamed.}
    \item{Typos in manual are corrected.}
  }
}

\section{Version 0.9.5 (2023-02-06)}{
  \itemize{
    \item{GLM, REG, aov1, aov2, aov3 functions return yhat and residuals with Resid=TRUE option.}
    \item{aov1, aov2, aov3 functions return beta coefficients with BETA=TRUE option.}
    \item{T3test is significantly improved for complex cases.}
    \item{RDmn handles more extreme values.}
    \item{LSM plotting order is ascending by default.}
  }
}

\section{Version 0.9.4 (2023-01-07)}{
  \itemize{
    \item{ORmn: Special cases of OR - 0, infinity or NaN - are checked in the input.}
  }
}

\section{Version 0.9.3 (2022-11-04)}{
  \itemize{
    \item{GLM and REG supports weighted regression with WEIGHTS option.}
    \item{REG displays tests using HC0, HC3, and White's first and second moment specification test.}
    \item{Code refactored to cope better with a singular or ill-conditioned matrix.}
    \item{Some unused arguments such as eps (in GLM, REG, aov1, aov2, and aov3) are removed.}
    \item{\code{e1}, \code{e2}, \code{e3} arguments are simplified to reduce unnecessary calculation.}
  }
}

\section{Version 0.9.2 (2022-10-05)}{
  \itemize{
    \item{corFisher function for correlation test by Fisher's Z transformation is added.}
    \item{RDmn, RRmn, ORmn use MN weights for point estimates also.}
    \item{RDmn1 finds confidence interval better.}
    \item{Some typo fixed.}
  }
}

\section{Version 0.9.1 (2022-09-05)}{
  \itemize{
    \item{RDinv, RRinv, ORinv, ORcmh for stratified or meta-analysis of 2x2 are added.}
    \item{GLM reports adjusted R-squared}
    \item{QuartileRange returns type 2 IQR by default.}
    \item{UNIV returns MAD and type 2 IQR too.}
    \item{Coll calculates like SAS, not like SPSS}
    \item{lfit help updated}
  }
}

\section{Version 0.9.0 (2022-08-05)}{
  \itemize{
    \item{RD, RR, OR, RDmn1, RRmn1, ORmn1, RDmn, RRmn, ORmn functions for 2x2 table with or without stratification are added.}
  }
}

\section{Version 0.8.1 (2022-05-10)}{
  \itemize{
    \item{Description corrected in satt help.}
  }
}

\section{Version 0.8.0 (2022-03-21)}{
  \itemize{
    \item{GLM function's arguments changed: SOLUTION to BETA, LSMEANS to EMEAN}
  }
}

\section{Version 0.7.0 (2022-03-14)}{
  \itemize{
    \item{ANOVA function is removed, because it is almost a duplicate of GLM function.}
    \item{GLM function's arguments changed.}
  }
}

\section{Version 0.6.7 (2022-03-12)}{
  \itemize{
    \item{Descriptive statistics functions such as N, Mean, SD, SEM, Min, Max, Median, UCL, LCL, trimmedMean, Range, QuartileRange, Skweness, SkewnessSE, Kurtosis, KurtosisSE now accept a non-numeric vector that can be converted to a numeric vector.}
    \item{A bug with no incercept model in \code{e3} is fixed.}
  }
}

\section{Version 0.6.6 (2022-02-19)}{
  \itemize{
    \item{A bug at "adj" option in \code{LSM} fixed, which also affects GLM.}
    \item{\code{PLOT} option related codes are revised.}
  }
}

\section{Version 0.6.5 (2022-01-20)}{
  \itemize{
    \item{\code{UNIV} added for descriptive statistics}
  }
}

\section{Version 0.6.4 (2021-11-25)}{
  \itemize{
    \item{\code{G2SWEEP} became faster for large sparse matrix}
  }
}

\section{Version 0.6.3 (2021-09-18)}{
  \itemize{
    \item{\code{REG} function warns in case of completely aliased model and data and increases epsilon tolerance to 1e-5.}
  }
}

\section{Version 0.6.2 (2021-08-19)}{
  \itemize{
    \item{\code{Diffogram} and \code{plotDiff} option 'Title' removed. Use 'main' option instead.}
  }
}

\section{Version 0.6.1 (2021-07-15)}{
  \itemize{
    \item{Descriptive statistic functions (Mean, SD, ...) return NA with inappropriate inputs.}
    \item{\code{geoMean} and \code{geoCV} added}
    \item{\code{reverse} option added to \code{PDIFF}}
  }
}

\section{Version 0.6.0 (2021-06-15)}{
  \itemize{
    \item{\code{Diffogram} added}
    \item{\code{PLOT} option added to \code{LSM} and \code{PDIFF}}
    \item{\code{LSM} supports grouping when a term is specified.}
    \item{\code{PDIFF} and \code{est} support common multiple comparison methods, such as "lsd"(no adj), "bon", "tukey", "scheffe", "duncan", and "dunnett".}
    \item{\code{ESTM} function added.}
    \item{\code{CONTR} function added.}
    \item{\code{SLICE} function added.}
    \item{\code{satt} changed the way of weights calculation. Now, one can give integer weights.}
    \item{\code{lsm} now properly supports models without intercept.}
    \item{\code{e2} now properly supports single factor models (Type II SS of Oneway ANOVA) without intercept.}
    \item{Checks inputs and errors during calculation more.}
  }
}

\section{Version 0.5.3 (2021-05-15)}{
  \itemize{
    \item{\code{aov2} and \code{e2} support no intercept model.}
    \item{\code{aov1}, \code{aov3}, \code{ANOVA}, \code{GLM} handle no intercept model properly.}
    \item{A column of estimability will be shown when some coefficients are not estimable in \code{REG} and \code{GLM}.}
  }
}

\section{Version 0.5.2 (2021-04-15)}{
  \itemize{
    \item{Some errata in help/manual are corrected.}
    \item{\code{CV} (coefficient of variation) in percentage added.}
    \item{\code{cSS} has now \code{m} arguments for non-zero mean hypothesis.}
    \item{A bug with no intercept model in \code{lr} fixed.}
    \item{Incomplete interaction term names of \code{PDIFF} fixed.}
  }
}

\section{Version 0.5.1 (2021-03-15)}{
  \itemize{
    \item{\code{LSM} does not show estimates which are not estimable, and supports more models}
    \item{\code{est} has more arguments to test estimability}
  }
}

\section{Version 0.5.0 (2021-02-15)}{
  \itemize{
    \item{\code{lr0} function for simple linear regressions with each independent variable}
  }
}

\section{Version 0.4.1 (2021-02-15)}{
  \itemize{
    \item{\code{lr} function for linear regression with g2 inverse}
    \item{\code{tsum} wrapper function for \code{tsum0}, \code{tsum1}, \code{tsum2}, \code{tsum3}}
    \item{\code{bt} function removed to remove dependency on \code{knitr}}
    \item{\code{NOINT} option in REG, ModelMatrix is removed. Use conventional R formula.}
    \item{Order of coefficients (parameters) changed to be more similar to SAS.}
    \item{\code{T3MS} revived.}
    \item{\code{EMS} has different arguments from T3MS.}
  }
}

\section{Version 0.4.0 (2021-01-29)}{
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

\section{Version 0.3.0 (2021-01-16)}{
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
    \item{TO DO: Function to add - e4, aov4}
    \item{TO DO: Items to add - R2, R2ADJ - v0.4.1 \code{lr} does this.}
  }
}
