w.sd =
  function (x, w)
  {
    sumw = sum(w)
    sumw.sq = sum(w)^2
    w.sq = sum(w^2)
    denom = sum(w * ((x - mean(x))^2))
    sqrt((sumw * denom) / (sumw.sq - w.sq))
  }
