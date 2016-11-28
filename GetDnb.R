
GetDnb <- function(ns, x, y,beta){

    uniq <- sort(unique(y))
    ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
    for(k in 1:length(uniq)){
      a <- colSums(x[y==uniq[k],])+beta
      b <- colSums(ns[y==uniq[k],])+beta
      ds[k,] <- a/b
    }
    return(ds)

}
