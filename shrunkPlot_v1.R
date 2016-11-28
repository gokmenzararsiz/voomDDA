shDiffPlot <- function(object, col, ...){
    nc <- ncol(object)
    nr <- nrow(object)
    
    xlims <- c(0.8, 2*nc + 1.2)
    ylims <- c(-0.2, nr + 0.2)
    
    par(mar = c(5.1,7,4.1,2.1))
    plot(c(max(xlims), min(xlims)), c(max(ylims), min(ylims)), type="n", axes = FALSE, ...)
    
    axis(1, at = seq(2,2*nc,2), labels = colnames(object))
    axis(2, at = 1:nr, labels = rownames(object), las = 2, cex.axis = 0.5)
    box()
    
    for (i in 1:nc){
        arrows(x0 = rep(2*i,nr), y0 = c(1:nr), x1 = 2*i + object[,i], y1 = c(1:nr) , length = 0, lwd=4, col = col[i])
        abline(v = 2*i, col = "gray70")
    }    
}