`intensityPlot` <-
function(figure, ## file name, such as "figure1.png"
                          data, ## log2 intensity file,
                          axisName,
                          mydir=getwd(),
                          figWidth= 540,
                          figHt = 360,
                          titleSize = 1.2,
                          axisSize = 0.7,
                          margin = c(6, 4, 3, 1), mycol){

    png(filename=paste(mydir, "/", figure, sep=""), width=figWidth, height=figHt)
    par(mar=margin)
    if (class(data) == "AffyBatch"){
        boxplot(data, names=axisName, main="Raw Intensity", 
            xlab="", ylab="", cex.main=titleSize, cex.axis=axisSize, cex.lab=1, 
            las=2, col=mycol)
    } else {
        boxplot(as.data.frame(data), names=axisName, main="Intensity Distribution", 
            xlab="", ylab="Log2 Intensity", cex.main=titleSize, cex.axis=axisSize, cex.lab=1, 
            las=2, col=mycol)
    }
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    dev.off()
}

