predict.voomNSC =
function (fit, newdata, threshold = NULL, prior = NULL) 
{
  if (is.null(prior)) prior = fit$prior
  if (is.null(threshold)) threshold = fit$opt.threshold
  vm = voomGSD(data.train = fit$counts, data.test = newdata, group = fit$conditions, norm.method = fit$normalization)
  x = vm$TestExp
  
  
# soft.shrink(...)
soft.shrink = 
  function (delta, threshold) 
  {
    dif = abs(delta) - threshold
    delta = sign(delta) * dif * (dif > 0)
    nonzero = sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
    attr(delta, "nonzero") = nonzero
    delta
  }

# diag.disc(...)
diag.disc =
  function (x, centroids, prior, weight) 
  {
    if (!missing(weight)) {
      posid = (weight > 0)
      if (any(posid)) {
        weight = sqrt(weight[posid])
        centroids = centroids[posid, , drop = FALSE] * weight
        x = x[posid, , drop = FALSE] * weight
      }
      else {
        mat = outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) = list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    dd = t(x) %*% centroids
    dd0 = drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - 
      log(prior)
    names(dd0) = NULL
    scale(dd, dd0, FALSE)
  }

# safe.exp(...)
safe.exp = function (x) 
{
  xx = sign(x) * pmin(abs(x), 500)
  return(exp(xx))
}

#soft.max(...)
softmax =
  function (x, gap = FALSE) 
  {
    d = dim(x)
    maxdist = x[, 1]
    pclass = rep(1, d[1])
    for (i in seq(2, d[2])) {
      l = x[, i] > maxdist
      pclass[l] = i
      maxdist[l] = x[l, i]
    }
    dd = dimnames(x)[[2]]
    if (gap) {
      x = abs(maxdist - x)
      x[cbind(seq(d[1]), pclass)] = drop(x %*% rep(1, d[2]))
      gaps = do.call("pmin", data.frame(x))
    }
    pclass = if (is.null(dd) || !length(dd)) 
      pclass
    else factor(pclass, levels = seq(d[2]), labels = dd)
    if (gap) 
      list(class = pclass, gaps = gaps)
    else pclass
  }

  sd = fit$weightedSD.pooled
  centroid.overall = fit$weightedMean
  centroids = fit$weightedMean.C
  se.scale = fit$se.scale
  delta = fit$delta
 
  delta.shrunk = soft.shrink(delta, threshold)
  delta.shrunk = t(t(delta.shrunk) * as.numeric(se.scale))  
  posid = drop(abs(delta.shrunk) %*% rep(1, length(prior))) > 0
  dd = diag.disc((x - fit$weightedMean)/(fit$weightedSD.pooled), delta.shrunk, prior, posid)
  softmax(dd)
}