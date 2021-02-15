BY = function(FUN, Formula, Data, By, ...)
{
  if (length(By) != 1) stop("Only 1 variable is supported!")
  if (!(By %in% colnames(Data))) stop(paste("Dataset does not have", By, "column!"))

  Levels = unique(Data[,By])
  nLevel = length(Levels)

  if (nLevel == 1) stop(paste(By, "has only 1 level!"))

  Result = list()
  for (i in 1:nLevel) {
    Result[[i]] = do.call(FUN, list(Formula, Data[Data[, By] == Levels[i],], ...))
  }

  names(Result) = Levels
  return(Result)
}
