w.mean =
  function (x, w) 
  {
    wx = w*x
    sum(wx) / sum(w)
  }