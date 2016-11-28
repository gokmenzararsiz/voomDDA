> diagDA
function (ls, cll, ts, pool = TRUE) 
{
  ls <- data.matrix(ls)
  n <- nrow(ls)
  p <- ncol(ls)
  cl0 <- as.integer(min(cll, na.rm = TRUE) - 1)
  cll <- as.integer(cll) - cl0
  inaC <- is.na(cll)
  clL <- cll[!inaC]
  K <- max(clL)
  if (K != length(unique(clL))) 
    stop(sQuote("cll"), " did not contain *consecutive* integers")
  nk <- integer(K)
  m <- v <- matrix(0, p, K)
  colVars <- function(x, means = colMeans(x, na.rm = na.rm), 
                      na.rm = FALSE) {
    x <- sweep(x, 2, means)
    colSums(x * x, na.rm = na.rm)/(nrow(x) - 1)
  }
  sum.na <- function(x) sum(x, na.rm = TRUE)
  for (k in 1:K) {
    which <- (cll == k)
    nk[k] <- sum.na(which)
    lsk <- ls[which, , drop = FALSE]
    m[, k] <- colMeans(lsk, na.rm = TRUE)
    if (nk[k] > 1) 
      v[, k] <- colVars(lsk, na.rm = TRUE, means = m[, 
                                                     k])
  }
  ts <- data.matrix(ts)
  if (p != ncol(ts)) 
    stop("test set matrix must have same columns as learning one")
  ts <- na.exclude(ts)
  nt <- nrow(ts)
  disc <- matrix(0, nt, K)
  if (pool) {
    vp <- rowSums(rep(nk - 1, each = p) * v)/(n - K)
    if (any(i0 <- vp == 0)) 
      vp[i0] <- 1e-07 * min(vp[!i0])
    ivp <- rep(1/vp, each = nt)
    for (k in 1:K) {
      y <- ts - rep(m[, k], each = nt)
      disc[, k] <- rowSums(y * y * ivp)
    }
  }
  else {
    if (FALSE) {
      for (k in 1:K) {
        ts <- ts - rep(m[, k], each = nt)
        disc[, k] <- rowSums((ts * ts)/rep(v[, k], each = nt)) + 
          sum(log(v[, k]))
      }
    }
    else {
      for (k in 1:K) {
        disc[, k] <- apply(ts, 1, function(z) sum((z - 
                                                     m[, k])^2/v[, k])) + sum.na(log(v[, k]))
      }
    }
  }
  pred <- cl0 + apply(disc, 1, which.min)
  if (inherits(attr(ts, "na.action"), "exclude")) 
    pred <- napredict(omit = attr(ts, "na.action"), pred)
  pred
}
<environment: namespace:sfsmisc>