predict.voomDDA = function(object, newdata){
  
  n = ncol(newdata)
  p = nrow(newdata)
  
  disc = matrix(0, n, object$nclass)  ## Discriminant scores for each class
  dimnames(disc) = list(colnames(newdata), object$classNames)
  vm = voomGSD(data.train = object$counts, data.test = newdata, group = object$conditions, norm.method = object$normalization)
  x = vm$TestExp

  x2 = t(x)

{
if (object$PooledVar) {
  vp = (object$weightedStats$weightedSD.pooled)^2
  if (any(i0 <- vp == 0)) 
    vp[i0] <- 1e-07 * min(vp[!i0])
    ivp <- rep(1/vp, each = n)
  for (k in 1:(object$nclass)) {
    y = x2 - rep(object$weightedStats$weightedMean.C[, k], each = n)
    disc[, k] = rowSums(y * y * ivp)
  }
}

else {
  if (FALSE) {
    for (k in 1:(object$nclass)) {
      x2 = x2 - rep(object$weightedStats$weightedMean.C[, k], each = n)
      vsd = (object$weightedStats$weightedSD.C)^2
      disc[, k] = rowSums((x2 * x2)/rep(vsd[, k], each = n)) + sum(log(vsd[, k]))
    }
  }
  else {
    vsd = (object$weightedStats$weightedSD.C)^2
    for (k in 1:(object$nclass)) {
      disc[, k] = apply(x2, 1, function(z) sum((z - object$weightedStats$weightedMean.C[, k])^2/vsd[, k])) + sum(log(vsd[, k]))
    }
  }
}
}

idx = apply(disc, 1, which.min)
pred = colnames(disc)[idx]
if (inherits(attr(x2, "na.action"), "exclude")) 
  pred = napredict(omit = attr(x2, "na.action"), pred)
pred
}