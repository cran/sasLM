## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(comment = NA)
options(width = 100)
library(sasLM)

## -------------------------------------------------------------------------------------------------
BEdata = af(BEdata, c("ADM", "SEQ", "PRD", "TRT", "SUBJ")) # columns as factors
formula1 = log(CMAX) ~ SEQ/SUBJ + PRD + TRT
GLM(formula1, BEdata)

## -------------------------------------------------------------------------------------------------
RanTest(formula1, BEdata, Random="SUBJ")

## -------------------------------------------------------------------------------------------------
ci0 = CIest(formula1, BEdata, "TRT", c(-1, 1), conf.level=0.90)
ci0
exp(ci0[, c("Estimate", "Lower CL", "Upper CL")]) # GMR and its 90% CI

## -------------------------------------------------------------------------------------------------
LSM(formula1, BEdata, "TRT", conf.level=0.90)

