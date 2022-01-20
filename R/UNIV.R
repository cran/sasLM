UNIV = function(x, conf.level = 0.95)
{
  x = as.numeric(x)
  if (!is.vector(x) | !is.numeric(x)) stop("A numeric vector is required!")
  nAll = length(x)
  nNA = sum(is.na(x))
  nFinite = sum(is.finite(x))
  if (nFinite == 0) stop("There is no finite number!")
  Res = c(nAll = nAll,
          nNA = nNA,
          nFinite = nFinite,
          Mean = Mean(x),
          Variance = var(x, na.rm=T),
          SD = SD(x),
          CV = CV(x),
          SEM = SEM(x),
          LowerConfLimit = LCL(x, conf.level=conf.level),
          UpperConfLimit = UCL(x, conf.level=conf.level),
          TrimmedMean = trimmedMean(x, Trim= 1 - conf.level),
          Min = Min(x),
          Q1 = quantile(x, 0.25, na.rm=T, names=F),
          Median = Median(x),
          Q3 = quantile(x, 0.75, na.rm=T, names=F),
          Max = Max(x),
          Range = Max(x) - Min(x),
          Skewness = Skewness(x),
          SkewnessSE = SkewnessSE(x),
          Kurtosis = Kurtosis(x),
          KurtosisSE = KurtosisSE(x))
  if (all(x > 0)) Res = c(Res, GeometricMean = geoMean(x), GeometricCV = geoCV(x))
  return(Res)
}

