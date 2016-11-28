pamr.train <-
function (data, n.threshold = 30, offset.percent = 50, remove.zeros = TRUE) 
{
  this.call <- match.call()
  if (!is.null(data$y)) {
    problem.type <- "class"
  }
 
  y <- as.factor(data$y)

  ytest <- NULL
  xtest <- NULL

  prior <- table(y)/length(y)

  junk <- nsc(data$x, y = y, offset.percent = offset.percent, n.threshold = n.threshold, prior = prior,
              remove.zeros = remove.zeros)
  
  junk$call <- this.call
  junk$problem.type <- problem.type
  class(junk) = "pamrtrained"
  junk
}
