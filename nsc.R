nsc <-
function (x, y, n.threshold = 30, offset.percent = 50, prior = NULL, remove.zeros = TRUE) 
{
  this.call <- match.call()
  Y <- model.matrix(~factor(y) - 1, data = list(y = y))
  
  xtest <- x
  ytest <- y
  
  
  n.class <- table(y)
  if (min(n.class) == 1) {
    stop(warning("Warning: a class contains only 1 sample"))
  }

  n <- sum(n.class)
  ntest <- ncol(xtest)
  K <- length(prior)
  p <- nrow(x)
  
  dimnames(Y) <- list(NULL, names(n.class))
  centroids <- scale(x %*% Y, FALSE, n.class)  ## WMEAN.G
  
  xdif <- x - centroids %*% t(Y)
  sd <- (xdif^2) %*% rep(1/(n - K), n)
  sd <- drop(sqrt(sd))  #WSD.POOLED
  offset <- quantile(sd, offset.percent/100)
  sd <- sd + offset

  centroid.overall <- drop(x %*% rep(1/n, n))   ## WMEAN
  
  se.scale <- sqrt(1/n.class - 1/n)  # mk
  
  delta <- (centroids - centroid.overall)/sd
  delta <- scale(delta, FALSE, se.scale)    ##dik
  
  threshold <- seq(0, max(abs(delta)), length = n.threshold)

  nonzero <- seq(n.threshold)
  errors <- threshold
  yhat <- as.list(seq(n.threshold))
  prob <- array(0, c(ntest, K, n.threshold))
  
  for (ii in 1:n.threshold) {
    cat(ii)
    delta.shrunk <- soft.shrink(delta, threshold[ii])
    #delta.shrunk <- scale(delta.shrunk, FALSE, 1/(se.scale))
    delta.shrunk <- t(t(delta.shrunk) * as.numeric(se.scale))
    
    nonzero[ii] <- attr(delta.shrunk, "nonzero")
    posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
    dd <- diag.disc((xtest - centroid.overall)/sd, delta.shrunk, 
                    prior, weight = posid)
    yhat[[ii]] <- softmax(dd)
    dd <- safe.exp(dd)
    prob[, , ii] <- dd/drop(dd %*% rep(1, K))
    if (!is.null(ytest)) {
      errors[ii] <- sum(yhat[[ii]] != ytest)
    }
  }
  thresh.names <- format(round(threshold, 3))
  names(yhat) <- thresh.names
  attr(yhat, "row.names") <- paste(seq(ntest))
  class(yhat) <- "data.frame"
  if (remove.zeros) 
    n.threshold <- match(0, nonzero, n.threshold)
  dimnames(prob) <- list(paste(seq(ntest)), names(n.class), 
                         thresh.names)
  object <- list(y = ytest, yhat = yhat, prob = prob[, , seq(n.threshold)], 
                 centroids = centroids, centroid.overall = centroid.overall, 
                 sd = sd, threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(n.threshold)], 
                 se.scale = se.scale, call = this.call, prior = prior, offset = offset)
  if (!is.null(ytest)) 
    object$errors <- errors
  #class(object) <- "nsc"
  object
}
