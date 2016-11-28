weighted.stats =
  function (x, w, conditions) 
  {
    n = ncol(x) #number of samples
    p = nrow(x) #number of genes
    nclass = length(unique(conditions)) #number of class
    if(is.factor(conditions)) {cNames = sort(levels(conditions))}
    if (is.numeric(conditions)) {cNames = as.character(sort(unique(conditions)))}
    
    WM = WS = wSum = se.scale = matrix(0, p, nclass)
    rownames(WS) = rownames(WM) = rownames(wSum) = rownames(se.scale) = rownames(x)
    colnames(WS) = colnames(WM) = colnames(wSum) = colnames(se.scale) = cNames
    
    c.ind = as.numeric(conditions)
    
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
      function (x, w, conditions) 
      {
        for (j in 1:nclass)
        {
          WM[,j] = w.mean00(x[,c.ind == j], w[,c.ind == j])
        }
        return(WM)
      }
    
    w.sd =
      function (x, w, conditions) 
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
        
        for (j in 1:nclass)
        {
          WS[,j] = w.sd00(x[,c.ind == j], w[,c.ind == j])
        }
        return(WS)
      }
    
    weightedMean = w.mean00(x, w) #Overall weighted mean
    weightedMean.C = w.mean(x, w, conditions) #Weighted means for each group 
    weightedSD.C = w.sd(x, w, conditions) #Weighted standard deviations for each group
    #weightedSD.pooled = weightedSD.C
    

    for (i in 1:nclass)
    {
      tmp = w[,which(c.ind == i)]
      rSum.tmp = rowSums(tmp)
      
      wSum[,i] = rSum.tmp
    }
    
    weightedSD.pooled = sqrt(rowSums((wSum-1) * (weightedSD.C^2)) / (rowSums(wSum) - nclass))
    se.scale = sqrt(1 / wSum + 1 / rowSums(wSum))
    
    s0 = median(weightedSD.pooled)
    
    delta = (weightedMean.C - weightedMean)/(se.scale*(weightedSD.pooled + s0))
        
    
    weightedSD.pooled = sqrt(rowSums(as.data.frame(weightedSD.pooled)) / (n - nclass))
    stats = list(n = n, p = p, nclass = nclass, se.scale = se.scale, weightedMean = weightedMean, weightedMean.C = weightedMean.C, weightedSD.C = weightedSD.C, weightedSD.pooled = weightedSD.pooled, delta = delta)
    
    return(stats)
  }