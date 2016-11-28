#Calculation of quantile normalization factors necessary for calcNormFactorsGSD function. Same for both train and test sets.
calcFactorQuantileGSD =
  function (data, lib.size, p = 0.75) 
  {
    y <- t(t(data)/lib.size)
    f <- apply(y, 2, function(x) quantile(x, p = p))
  }

#Calculation of weighted normalization factors necessary for calcNormFactorsGSD function.
calcFactorWeightedGSD =
  function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3, 
            sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10) 
  {
    if (all(obs == ref)) 
      return(1)
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    if (is.null(libsize.obs)) 
      nO <- sum(obs)
    else nO <- libsize.obs
    if (is.null(libsize.ref)) 
      nR <- sum(ref)
    else nR <- libsize.ref
    logR <- log2((obs/nO)/(ref/nR))
    absE <- (log2(obs/nO) + log2(ref/nR))/2
    v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    n <- sum(fin)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)
    if (doWeighting) 
      2^(sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],na.rm = TRUE))
    else 2^(mean(logR[keep], na.rm = TRUE))
  }

#Calculation of deseq normalization factors necessary for calcNormFactorsGSD function.
calcFactorRLEGSD =
  function (data.train, data.test, lib.size, lib.size.test) 
  {
    gm <- exp(rowMeans(log(data.train)))
    f = apply(data.train, 2, function(u) median((u/gm)[gm > 0]))
    f.test = apply(data.test, 2, function(u) median((u/gm)[gm > 0]))
    f = f / lib.size
    f.test = f.test / lib.size.test
    deseqsizefactors = list(f,f.test)
    return(deseqsizefactors)
  }

#Calculation of normalization factors using calcNormFactorsGSD function.
calcNormFactorsGSD =
  function (data.train, data.test, lib.size = NULL, method = c("TMM", "deseq", "none"), refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, 
            doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, ...) 
  {
    x <- as.matrix(data.train)
    xtest <- as.matrix(data.test)
    if (any(is.na(x)||is.na(xtest)))
      stop("NAs not permitted")
    if (is.null(lib.size)) 
      lib.size <- colSums(x)
      
    lib.size.test <- colSums(xtest)
    method <- match.arg(method)
    allzero <- rowSums(x > 0) == 0
    if (any(allzero)) 
      x <- x[!allzero, , drop = FALSE]
      xtest <- xtest[!allzero, , drop = FALSE]
    if (nrow(x) == 0 || ncol(x) == 1) 
      method = "none"
    
    if(method == "TMM"){
      f75 <- calcFactorQuantileGSD(data = x, lib.size = lib.size, p = 0.75)
      f75.test <- calcFactorQuantileGSD(data = xtest, lib.size = lib.size.test, p = 0.75)
      
      refColumn <- which.min(abs(f75 - mean(f75)))
      f <- rep(NA, ncol(x))
      f.test <- rep(NA, ncol(xtest))
      for (i in 1:ncol(x)) f[i] <- calcFactorWeightedGSD(obs = x[,i], ref = x[, refColumn], libsize.obs = lib.size[i], 
                                                         libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
                                                         sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
      for (i in 1:ncol(xtest)) f.test[i] <- calcFactorWeightedGSD(obs = xtest[,i], ref = x[, refColumn], libsize.obs = lib.size.test[i], 
                                                                  libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
                                                                  sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
      normf = list(f,f.test)
    } 
    else if(method == "deseq"){
      normf = calcFactorRLEGSD(data.train = x, data.test = xtest, lib.size = lib.size, lib.size.test = lib.size.test)#/lib.size 
    }
    else {
      normf = list(rep(1, ncol(x)), rep(1, ncol(xtest)))
    }
    
    names(normf) = c("train", "test")
    
    f = as.numeric(normf[[1]]) / (exp(mean(log(normf[[1]]))))
    f.test = as.numeric(normf[[2]]) / (exp(mean(log(normf[[1]]))))
    normf2 = list(f, f.test, lib.size, lib.size.test)
    names(normf2) = c("TrainNormFactor","TestNormFactor","TrainLibSize","TestLibSize")
    return(normf2)
  }

#traindatayÄ±, testdatayÄ± ve trainclassÄ± alacak, method = TMM, RLE, none olacak.
#train ve test icin voom gibi expression ve weights dÃ¶ndÃ¼recek.
voomGSD = 
  function(data.train, data.test, group, norm.method = c("TMM", "deseq", "none"), design = NULL, lib.size = NULL, span = 0.5)
  {
    out <- list()
    NormFactors = calcNormFactorsGSD(data.train = data.train, data.test = data.test, method = norm.method)
    TrainNormFactor = NormFactors$TrainNormFactor
    TestNormFactor = NormFactors$TestNormFactor
    TrainLibSize = NormFactors$TrainLibSize
    TestLibSize = NormFactors$TestLibSize
    lib.size.tr = TrainNormFactor * TrainLibSize
    lib.size.ts = TestNormFactor * TestLibSize
    
    design.tr = model.matrix(~group)
    rownames(design.tr) = colnames(data.train)
    
    design.ts <- matrix(1, ncol(data.test), 1)
    rownames(design.ts) <- colnames(data.test)
    colnames(design.ts) <- "GrandMean"
    
    y.tr <- t(log2(t(data.train + 0.5)/(lib.size.tr + 1) * 1e+06))
    y.ts <- t(log2(t(data.test + 0.5)/(lib.size.ts + 1) * 1e+06))
    fit.tr <- lmFit(y.tr, design.tr)
    fit.ts <- lmFit(y.ts, design.ts)

    
    if (is.null(fit.tr$Amean)) 
      fit$Amean <- rowMeans(y.tr, na.rm = TRUE)
    
    fit.ts$Amean = fit.tr$Amean
    fit.ts$sigma = fit.tr$sigma
    fit.ts$coefficients = fit.tr$coefficients[,1]
    
    sx <- fit.tr$Amean + mean(log2(lib.size.tr + 1)) - log2(1e+06)
    sy <- sqrt(fit.tr$sigma)
    l <- lowess(sx, sy, f = span)
    f <- approxfun(l, rule = 2)

    fitted.values.tr <- fit.tr$coefficients %*% t(fit.tr$design)
    fitted.values.ts <- fit.ts$coefficients %*% t(fit.ts$design)
    fitted.cpm.tr <- 2^fitted.values.tr
    fitted.cpm.ts <- 2^fitted.values.ts
    fitted.count.tr <- 1e-06 * t(t(fitted.cpm.tr) * (lib.size.tr + 1))
    fitted.count.ts <- 1e-06 * t(t(fitted.cpm.ts) * (lib.size.ts + 1))
    fitted.logcount.tr <- log2(fitted.count.tr)
    fitted.logcount.ts <- log2(fitted.count.ts)
    w.tr <- 1/f(fitted.logcount.tr)^4
    w.ts <- 1/f(fitted.logcount.ts)^4
    dim(w.tr) <- dim(fitted.logcount.tr)
    dim(w.ts) <- dim(fitted.logcount.ts)
    dimnames(w.tr) = dimnames(y.tr)
    dimnames(w.ts) = dimnames(y.ts)
    out$TrainExp <- y.tr
    out$TestExp <- y.ts
    out$TrainWeights <- w.tr
    out$TestWeights <- w.ts
    new("EList", out)
  }
