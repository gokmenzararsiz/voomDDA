weighted.stats =
  function (x, w, c) 
  {
    n = ncol(x) #number of samples
    p = nrow(x) #number of genes
    k = length(unique(c)) #number of class
    c = as.integer(c)
    WM = WS = matrix(0, p, k)
    rownames(WS) = rownames(x)
    colnames(WS) = unique(c)
    c0 <- as.integer(min(c, na.rm = TRUE) - 1)
    c <- as.integer(c) - c0
    mk = NULL
    
    w.mean00 =
      function (x, w) 
      {
        wm = NULL
        
        for (i in 1:p)
        {
          wm0 = sum(w[i,]*x[i,]) / sum(w[i,])
          wm = c(wm, wm0)
        }
        return(wm)
      }
    
    w.mean =
      function (x, w, c) 
      {
        for (j in 1:k)
        {
          WM[,j] = w.mean00(x[,c == j], w[,c == j])
        }
        return(WM)
      }
    
    w.sd =
      function (x, w, c) 
      {
        w.sd00 =
          function (x, w) 
          {
            ws = NULL
            
            w.sd0 =
              function (x, w)
              {
                sumw = sum(w)
                sumw.sq = sum(w)^2
                w.sq = sum(w^2)
                denom = sum(w * ((x - mean(x))^2))
                sqrt((sumw * denom) / (sumw.sq - w.sq))
              }
            
            for (i in 1:p)
            {
              ws0 = w.sd0(x[i,], w[i,])
              ws = c(ws, ws0)
            }
            
            return(ws)
          }
        
        for (j in 1:k)
        {
          WS[,j] = w.sd00(x[,c == j], w[,c == j])
        }
        return(WS)
      }
    
    WMEAN = w.mean00(x, w) #Overall weighted mean
    delta = WMEAN.G = w.mean(x, w, c) #Weighted means for each group 
    WSD.G = w.sd(x, w, c) #Weighted standard deviations for each group
    WSD.POOLED = WSD.G
    
    for (i in 1:k)
    {
      WSD.POOLED[,i] = (table(c)[i]-1) * (WSD.POOLED[,i]^2)
      mk[i] = sqrt((1 / table(c)[i]) + (1 / n))
    }
    
    WSD.POOLED = sqrt(rowSums(as.data.frame(WSD.POOLED)) / (n - k))
    s0 = median(WSD.POOLED)
    
    for (i in 1:k)
    {
      delta[,i] = (delta[,i] - WMEAN) / (mk[i]*(WSD.POOLED + s0))
    }
    
    rownames(WMEAN) = rownames(WMEAN.G) = rownames(delta) = rownames(WSD.G) = rownames(WSD.POOLED) = rownames(x)
    colnames(WMEAN.G) = colnames(delta) = colnames(WSD.G) = unique(c)
    
    stats = list(n = n, p = p, k = k, mk = mk, s0 = s0, delta = delta, WMEAN = WMEAN, WMEAN.G = WMEAN.G, WSD.G = WSD.G, WSD.POOLED = WSD.POOLED)
    
    return(stats)
  }