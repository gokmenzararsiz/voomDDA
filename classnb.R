classnb <- function (x, y, xte = NULL, beta = 1,  type = c("mle", 
    "deseq", "quantile"), prior = NULL) 
{

    if (is.null(prior)) 
        prior <- rep(1/length(unique(y)), length(unique(y)))

    null.out <- NullModel(x, type = type)
    ns <- null.out$n
    nste <- NullModelTest(null.out, x, xte, type = type)$nste

    uniq <- sort(unique(y))

        ds <- GetDnb(ns, x, y, beta)

        #disperhat=rep(truephi,ncol(nste))
        phihat <- as.numeric(disperhat)
        discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

        for (k in 1:length(uniq)) {

            for(l in 1:nrow(xte))   {

                 dstar = ds[k,]
                 part2=1+nste[l,]*dstar*phihat 
                 part1=dstar/part2 

                 discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])

             }
         }
        save <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, 
            ytehat = uniq[apply(discriminant, 1, which.max)], 
            x = x, y = y, xte = xte, 
            type = type)

        return(save)

}




