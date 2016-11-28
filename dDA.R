
> dDA
function (x, cll, pool = TRUE) 
{
  x <- data.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
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
    lsk <- x[which, , drop = FALSE]
    m[, k] <- colMeans(lsk, na.rm = TRUE)
    if (nk[k] > 1) 
      v[, k] <- colVars(lsk, na.rm = TRUE, means = m[, 
                                                     k])
  }
  structure(list(call = match.call(), cl0 = cl0, n = n, p = p, 
                 K = K, means = m, vars = v, nk = nk, pool = pool), class = "dDA")
}
