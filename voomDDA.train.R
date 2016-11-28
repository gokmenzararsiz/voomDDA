voomDDA.train = function (counts, conditions, normalization = c("TMM", "deseq", "none"), pooled.var = TRUE) 
{
  normalization = match.arg(normalization)
  counts = data.matrix(counts)
  p = nrow(counts)
  n = ncol(counts)
  
  nclass = length(unique(conditions))
  n.class = table(conditions)
  
  if (min(n.class) == 1) {
    stop(warning("Warning: a class contains only 1 sample"))
  }
 
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
  
  if(is.factor(conditions)) classNames = levels(conditions)

  wStats = weighted.stats(x, w, conditions)
  
  results = structure(list(counts = counts, conditions = factor(conditions), weightedStats = wStats, normalization = normalization,
                           nclass = nclass, classNames = classNames, PooledVar = pooled.var), class = "voomDDA")
  
  return(results)
}