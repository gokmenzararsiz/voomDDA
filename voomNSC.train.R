voomNSC.train =
  function (counts, conditions, n.threshold = 30, offset.percent = 50, remove.zeros = TRUE,
            normalization = c("TMM", "deseq", "none")) 
  {
    normalization = match.arg(normalization)
    this.call = match.call()
    
#### voom transformation steps
    {
      if (normalization == "TMM"){
        design = model.matrix(~conditions)
        rownames(design) = colnames(counts)
        dge = DGEList(counts = counts)
        dge = calcNormFactors(dge, method = "TMM")
        vm = voom(dge, design, plot=F)
        x = vm $ E #pxn dim gene expression matrix
        w = vm $ weights #pxn dim weight matrix
        dimnames(w) = dimnames(x)
      } 
      else if (normalization == "deseq"){
        design = model.matrix(~conditions)
        rownames(design) = colnames(counts)
        dge = DGEList(counts = counts)
        dge = calcNormFactors(dge, method = "RLE")
        vm = voom(dge, design, plot=F)
        x = vm $ E #pxn dim gene expression matrix
        w = vm $ weights #pxn dim weight matrix
        dimnames(w) = dimnames(x)
      }
      else {
        design = model.matrix(~conditions)
        rownames(design) = colnames(counts)
        vm = voom(counts, design, plot=F)
        x = vm $ E #pxn dim gene expression matrix
        w = vm $ weights #pxn dim weight matrix
        dimnames(w) = dimnames(x)
      }
    }

    y = as.factor(conditions)
    n.class = table(y)
    prior = n.class/length(y)
    
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

## wnsc(....)
wnsc =
  function (x, y, n.threshold = 30, offset.percent = 50, prior = NULL, remove.zeros = TRUE, weights = NULL) 
  {
    selected.genes = selected.genesIndex = list()
    this.call = match.call()
    Y = model.matrix(~factor(y) - 1, data = list(y = y))
    
    xtest = x
    ytest = y
    
    wStats = weighted.stats(x, weights, y)
    
if (min(n.class) == 1) {
      stop(warning("Warning: a class contains only 1 sample"))
    }
    
    n = sum(n.class)
    ntest = ncol(xtest)
    K = length(prior)
    p = nrow(x)
    dimnames(Y) = list(NULL, names(n.class))
    centroids = wStats$weightedMean.C
    sd = wStats$weightedSD.pooled
    offset = quantile(sd, offset.percent/100)
    sd = sd + offset
    centroid.overall = wStats$weightedMean
    se.scale = wStats$se.scale
    delta = wStats$delta
    threshold = seq(0, max(abs(delta)), length = n.threshold)
    nonzero = seq(n.threshold)
    errors = threshold
    pred.conditions = as.list(seq(n.threshold))
    prob = array(0, c(ntest, K, n.threshold))
    
    dshrunkAll <- list()
    for (ii in 1:n.threshold) {
      
      delta.shrunk = soft.shrink(delta, threshold[ii])
      delta.shrunk = t(t(delta.shrunk) * as.numeric(se.scale))
      dshrunkAll[[ii]] <- delta.shrunk
      
      nonzero[ii] = attr(delta.shrunk, "nonzero")
      posid = drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      selected.genes[[ii]] = rownames(wStats$weightedMean.C)[posid]
      selected.genesIndex[[ii]] = which(posid == TRUE)
      dd = diag.disc((xtest - centroid.overall)/sd, delta.shrunk, 
                      prior, weight = posid)
      pred.conditions[[ii]] = softmax(dd)
      dd = safe.exp(dd)
      prob[, , ii] = dd/drop(dd %*% rep(1, K))
      if (!is.null(ytest)) {
        errors[ii] = sum(pred.conditions[[ii]] != ytest)
      }
    }
    
    thresh.names = format(round(threshold, 3))
    names(pred.conditions) = names(selected.genes) = thresh.names
    attr(pred.conditions, "row.names") = paste(seq(ntest))
    class(pred.conditions) = "data.frame"
    
    if (remove.zeros) 
      n.threshold = match(0, nonzero, n.threshold)
    
    dimnames(prob) = list(paste(seq(ntest)), names(n.class), 
                           thresh.names)
    
    object = list(counts = counts, conditions = factor(conditions), pred.conditions = pred.conditions, prob = prob[, , seq(n.threshold)], 
                  weightedMean.C = centroids, weightedMean = centroid.overall, delta = delta, normalization = normalization,
                  weightedSD.pooled = sd, threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(n.threshold)], 
                   se.scale = se.scale, call = this.call, prior = prior, offset = offset, SelectedGenes = selected.genes, 
                   SelectedGenesIndex = selected.genesIndex)
    
    object$errors = errors
    
    opt.threshold = max(object$threshold[which(object$errors == min(object$errors))])
    
    model.res = as.data.frame(cbind(object$threshold, object$nonzero, object$errors))
    colnames(model.res) = c("threshold", "nonzero", "errors")
    object$modelRes = model.res
    
    object$delta.shrunk <- dshrunkAll[[which(object$threshold == opt.threshold)]]
    #min.error = model.res[model.res$errors == min(model.res$errors),]
    #min.error.gene = min.error[min.error$nonzero == min(min.error$nonzero),]
    #opt.threshold = min.error.gene$threshold
    
    object$opt.threshold = opt.threshold
    object
  }

    junk = wnsc(x, y = conditions, offset.percent = offset.percent, n.threshold = n.threshold, prior = prior,
                remove.zeros = remove.zeros, weights = w)
    
    junk$call = this.call
    class(junk) = "pamrtrained"
    junk
  }