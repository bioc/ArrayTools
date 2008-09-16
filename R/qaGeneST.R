`qaGeneST` <-
function(object, parameters, QC, mydir=getwd(), 
    outputFile = "QA.html"){
    
    require("R2HTML")
    targetsFile <- pData(object)
    expr <- exprs(object)

    num <- match(parameters, colnames(targetsFile))
    if (length(num) < 1 | sum(is.na(num)) > 0) 
        stop("Please check your parameter names, which must match the 
        \n  colnames of phenotype files.")
    targetsFile <- targetsFile[, num, drop=F]
    nchip <- ncol(expr)    
   
    targetSort<-targetsFile[do.call("order", targetsFile[,num, drop=F]),, drop=F]
   
    temp <- as.vector(table(targetSort))
    targetSort$SampleID<- unlist(lapply(temp[temp !=0], seq))
    targetTag <- apply(sapply(targetSort, as.character), 1, paste, collapse=".")
    targetSort$SampleID<- targetTag   ## index column
    logdata4QC <- expr[,match(rownames(targetSort), colnames(expr))]

    ##HTMLStart(outdir =mydir, filename="QA")
    cur.dir <- mydir
    result.dir <- paste(cur.dir, "Result", sep="/")
    dir.create(path=result.dir, showWarnings = FALSE)
    setwd(result.dir)
    QA.dir <- paste(result.dir, "Quality_Control_Assessment", sep="/")
    dir.create(path=QA.dir, showWarnings = FALSE)
    setwd(QA.dir)

    
    plotQC <- function(figure, QC, axisName, index, title, mycol, mydir=getwd(),
        margin=c(6, 4, 3, 1), figWidth= 540, figHt = 360, titleSize = 1.2, axisSize = 0.7, 
        legd=TRUE, logvalue=FALSE) { 
     
        nchip <- nrow(QC)
    
 
        if (logvalue) {
            QC1 <- log2(QC[,index])
        } else {
            QC1 <- QC[,index]
        }

        png(filename=paste(mydir, "/", figure, sep=""), width=figWidth, height=figHt)

        ymin <- min(QC1) / 1.25
        ymax <- ifelse (legd==TRUE, max(QC1) * 1.5, max(QC1) * 1.25)
        par(xaxt="n",mar=margin , las=2, cex.main=titleSize, cex.lab=1, cex.axis=axisSize)    
        if (length(index) == 1) {
      plot(as.numeric(QC1), ylim=c(ymin, ymax), pch=1, 
                col=mycol[1], type="b", ylab="", xlab="", cex=1.5, lty=3, lwd=3) 
            par(xaxt="s")
            axis(side=1, at=1:nchip, labels=axisName)  
        }  else {
            for (i in 1:length(index)) {
    plot(as.numeric(QC1[,i]), ylim=c(ymin, ymax), pch=i, 
                    col=rainbow(length(index))[i], type="b", ylab="", xlab="", cex=1.5, 
                    lty=i, lwd=2)
    par(new=T)
      }
            par(xaxt="s")
            axis(side=1, at=1:nchip, labels=axisName)  
            if (legd== TRUE){
                legend(x=1, y=ymax*0.98, legend=colnames(QC1), lty=1:length(index), 
                    pch=1:length(index), col=rainbow(length(index)), 
ncol=2, cex=0.8)
            }
        }     
        title(main=toupper(title), xlab="", ylab="Value")
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        dev.off()
    }

    HTMLStart(filename=outputFile)
    ##HTMLStart(filename=paste(mydir, outputFile, sep="\\"))

    titl <- as.title("<p align=center>Quality Assessment</p>")
    HTML(titl, file = outputFile , append = FALSE)
    HTML(targetSort, file = outputFile, Border=1, innerBorder=1)
    ##HTML(paste("The samples are grouped by:", paste(parameters, 
    ##    collapse = " and ")), file = outputFile )
    HTML("<li> Figure 1 -<a href= #fig1> Intensity Distribution </a>", 
        file = outputFile )
    HTML("<li> Figure 2 -<a href= #fig2> Mean Signal </a>", 
        file = outputFile )
    HTML("<li> Figure 3 -<a href= #fig3> BAC Spike </a>", 
        file = outputFile )
    HTML("<li> Figure 4 -<a href= #fig4> polya Spike </a>", 
        file = outputFile )
    HTML("<li> Figure 5 -<a href= #fig5> Pos Vs Neg Auc </a>", 
        file = outputFile )
    HTML("<li> Figure 6 -<a href= #fig6> Mad Residual Signal </a>", 
        file = outputFile )
    HTML("<li> Figure 7 -<a href= #fig7> RLE MEAN </a>", 
        file = outputFile )
    HTML("<li> Figure 8 -<a href= #fig8> Hierarchical Clustering of Samples </a>", 
        file = outputFile )
 
    figure1 <- "figure1.png"
    ##browser()
    HTML("<hr><a name= \"fig1\"></a>", file = outputFile)
    intensityPlot(figure1, logdata4QC, targetTag, mycol=rainbow(nchip), mydir=QA.dir)    
    HTMLInsertGraph(figure1, file = outputFile, Caption = "Figure 1: Intensity 
        Distribution", GraphBorder = 1, Align = "center", append = TRUE, 
        WidthHTML=540, HeightHTML=360)

    figure2 <- "figure2.png"
    rle_index <- grep("rle_mean", colnames(QC))
    mad_index <- grep("mad_residual_mean", colnames(QC))
    mean_item <- grep("_mean", colnames(QC))
    exclude_item <- c(which(colnames(QC)=="bac_spike_mean"), 
        which(colnames(QC)=="polya_spike_mean"))
    mean_index <- setdiff(mean_item, c(rle_index, mad_index, exclude_item))
    HTML("<hr><a name= \"fig2\"></a>", file = outputFile)
    plotQC(figure2, QC, targetTag , mean_index, "Mean Signal", mydir=QA.dir, logvalue=TRUE)
    HTMLInsertGraph(figure2, file = outputFile, Caption = "Figure 2: Mean Signal", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure3 <- "figure3.png"
    bac_5_at <- grep("5_at$", colnames(QC))
    polya_spike <- grep("_st$", colnames(QC))
    HTML("<hr><a name= \"fig3\"></a>", file = outputFile, logvalue=TRUE)
    plotQC(figure3, QC, targetTag , bac_5_at, "BAC Spike", mydir=QA.dir, logvalue=TRUE)
    HTMLInsertGraph(figure3, file = outputFile, Caption = "Figure 3: BAC Spike", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure4 <- "figure4.png"
    HTML("<hr><a name= \"fig4\"></a>", file = outputFile)
    plotQC(figure4, QC, targetTag , polya_spike, "polya Spike", mydir=QA.dir, logvalue=TRUE)
    HTMLInsertGraph(figure4, file = outputFile, Caption = "Figure 4: polya Spike", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure5 <- "figure5.png"
    HTML("<hr><a name= \"fig5\"></a>", file = outputFile)
    plotQC(figure5, QC, targetTag , grep("pos_vs_neg_auc", colnames(QC)), 
        "pos vs neg auc", mycol="green", mydir=QA.dir, legd=FALSE)
    HTMLInsertGraph(figure5, file = outputFile, Caption = "Figure 5: Pos Vs Neg Auc", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)
    
    figure6 <- "figure6.png"
    HTML("<hr><a name= \"fig6\"></a>", file = outputFile)
    plotQC(figure6, QC, targetTag , mad_index, "mad residual mean", mydir=QA.dir)
    HTMLInsertGraph(figure6, file = outputFile, Caption = "Figure 6: Mad Residual Signal", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure7 <- "figure7.png"
    HTML("<hr><a name= \"fig7\"></a>", file = outputFile)
    plotQC(figure7, QC, targetTag , rle_index, "rle mean", mydir=QA.dir)
    HTMLInsertGraph(figure7, file = outputFile, Caption = "Figure 7: RLE MEAN", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure8 <- "figure8.png"
    HTML("<hr><a name= \"fig8\"></a>", file = outputFile)
    plotCluster (figure8, logdata4QC, targetTag , "Hierarchical Clustering of Samples",
        mydir=QA.dir)
    HTMLInsertGraph(figure8, file = outputFile, Caption = "Figure 8: Hierarchical 
        Clustering of Samples", GraphBorder = 1, Align = "center", append = TRUE, 
        WidthHTML=540, HeightHTML=360)
    HTMLStop()
    setwd(cur.dir)
}

