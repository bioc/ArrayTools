`plotCluster` <-
function(figure, data, axisName, title, mydir=getwd(),
    figWidth= 540, figHt = 360, margin=c(2, 4, 3, 1)){
    corr <- cor(data)
    hc <- hclust(as.dist(1-corr), method="average")
    png(filename=paste(mydir, "/", figure, sep=""), width=figWidth, height=figHt)
    par(mar=margin)
    plot(hc, main = title, xlab = " ", sub = " ", labels = axisName, lwd=2)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    dev.off()
}

