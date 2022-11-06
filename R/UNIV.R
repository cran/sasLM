UNIV = function(y, conf.level = 0.95)
{
  y = as.numeric(y)
  if (!is.vector(y) | !is.numeric(y)) stop("A numeric vector is required!")
  nAll = length(y)
  nNA = sum(is.na(y))
  nFinite = sum(is.finite(y))
  if (nFinite == 0) stop("There is no finite number!")
  Res = c(nAll = nAll,
          nNA = nNA,
          nFinite = nFinite,
          Mean = Mean(y),
          Variance = var(y, na.rm=T),
          SD = SD(y),
          CV = CV(y),
          SEM = SEM(y),
          LowerConfLimit = LCL(y, conf.level=conf.level),
          UpperConfLimit = UCL(y, conf.level=conf.level),
          TrimmedMean = trimmedMean(y, Trim= 1 - conf.level),
          Min = Min(y),
          Q1 = quantile(y, 0.25, na.rm=T, names=F, type=6),
          Median = Median(y),
          Q3 = quantile(y, 0.75, na.rm=T, names=F, type=6),
          Max = Max(y),
          Range = Max(y) - Min(y),
          IQR = IQR(y, na.rm=T, type=6), # SAS default
          MAD = mad(y, na.rm=T),
          Skewness = Skewness(y),
          SkewnessSE = SkewnessSE(y),
          Kurtosis = Kurtosis(y),
          KurtosisSE = KurtosisSE(y))
  if (all(y > 0)) Res = c(Res, GeometricMean = geoMean(y), GeometricCV = geoCV(y))
  return(Res)
}

