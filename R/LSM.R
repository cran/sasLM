LSM = function(Formula, Data, conf.level=0.95)
{
  x = ModelMatrix(Formula, Data)
  y = model.frame(Formula, Data)[,1]
  rx = lfit(x, y)

  Res = lsm0(x, rx, Formula, Data, conf.level = conf.level)
  return(Res)
}

