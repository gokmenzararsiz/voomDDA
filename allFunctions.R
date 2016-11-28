
fitModel <- function(iter = 5, ...){
    error = TRUE
    i = 1
    while ((error) && (i < 5)){
        model <- try(train(...))
        error <- class(model) == "try-error"
        i = i + 1
    }
    if (i == iter) {
        model = NA
        warning("Algorithm did not converge. Returning NAs")
    }
    return(model)
}

predictModel <- function(model, newdata){
    if (class(model) == "train"){
        predicted = predict.train(model, newdata = newdata)    
        ts.acc = sum(test.cond == as.numeric(predicted)) / length(test.cond)
    }
    if (class(model) != "train"){ 
        ts.acc = NA
    }
    return(ts.acc)
}


foldIndex <- function(data, nSamp = NULL, nFolds = 5, repeats = 2){
    if(is.null(nSamp)) n = nrow(data)
    else n = nSamp
    
    indIn <- indOut <- list()
    
    for (j in 1:repeats){
        tmpIn = createFolds(1:n, k = nFolds, list = TRUE, returnTrain = TRUE)
        tmpOut = lapply(tmpIn, function(x)c(1:n)[-x])
        
        indIn = c(indIn, tmpIn)
        indOut = c(indOut, tmpOut)
    }
    
    nms = paste(rep(paste("Fold", 1:nFolds, sep = ""), repeats), 
                rep(paste(".Rep", 1:repeats, sep = ""), c(rep(nFolds, repeats))), sep = "")
    names(indIn) <- names(indOut) <- nms
    return(list(indexIn = indIn, indexOut = indOut))
}

