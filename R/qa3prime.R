`qa3prime` <-
function (object, parameters, outputFile = "QA.html", mydir=getwd()) 
{
    ##browser()
    require("R2HTML")
    require("affyPLM")
    require("simpleaffy")

    targetsFile <- pData(object)
    num <- match(parameters, colnames(targetsFile))
    if (length(num) < 1 | sum(is.na(num)) > 0) 
        stop("Please check your parameter names, which must match the 
        \n  colnames of phenotype files.")
    targetsFile <- targetsFile[, num, drop=F]
    nchip <- nrow(targetsFile )    
   
    targetSort<-targetsFile[do.call("order", targetsFile[,num, drop=F]),, drop=F]
    temp <- as.vector(table(targetSort))
    targetSort$SampleID<- unlist(lapply(temp[temp !=0], seq))
    targetTag <- apply(sapply(targetSort, as.character), 1, paste, collapse=".")
    targetSort$SampleID<- targetTag   ## index column

    object <- object[, match(rownames(targetSort), rownames(pData(object)))]
    genome <- cdfName(object)
    testchip <- grep("^Test3", genome)
    if (length(testchip) > 0) {
        genome <- "test3"
        object@cdfName <- "test3"
    }
    genome <- cleancdfname(cdfName(object))
    setQCEnvironment(genome)
    included <- qc.get.alpha1()

    ## Set directory to save files
    cur.dir <- mydir
    result.dir <- paste(cur.dir, "Result", sep="/")
    dir.create(path=result.dir, showWarnings = FALSE)
    setwd(result.dir)
    QA.dir <- paste(result.dir, "Quality_Control_Assessment", sep="/")
    dir.create(path=QA.dir, showWarnings = FALSE)
    setwd(QA.dir)
    
    emptyPlot <-function(figureName){
        png(filename = figureName)
        plot(1, 1, col = "white", yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
        text(1, 1, label = "Your Chip is not supported by simplyaffy", col = "red")
        dev.off()
    }

    lineGraph <- function(figure, QC, axisName, title, mydir=getwd(),
        margin=c(6, 4, 3, 1), figWidth= 540, figHt = 360, titleSize = 1.2, axisSize = 0.7, 
        legd=TRUE, logvalue=FALSE) {   
        nchip <- nrow(QC) 
        png(filename=paste(mydir, "/", figure, sep=""), width=figWidth, height=figHt)
        ymin <- min(QC) / 1.25
        ymax <- ifelse (legd==TRUE, max(QC) * 1.5, max(QC) * 1.25)
        par(xaxt="n",mar=margin , las=2, cex.main=titleSize, cex.lab=1, cex.axis=axisSize)    
        if (ncol(QC) == 1) {
      plot(as.numeric(QC[,1]), ylim=c(ymin, ymax), pch=1, 
                col=rainbow(1), type="b", ylab="", xlab="", cex=1.5, lty=3, lwd=3) 
            par(xaxt="s")
            axis(side=1, at=1:nchip, labels=axisName)  
        }  else {
            for (i in 1:ncol(QC)) {
             plot(as.numeric(QC[,i]), ylim=c(ymin, ymax), pch=i, 
                    col=rainbow(ncol(QC))[i], type="b", ylab="", xlab="", cex=1.5, 
                    lty=i, lwd=2)
       par(new=T)
           }
            par(xaxt="s")
            axis(side=1, at=1:nchip, labels=axisName)  
            if (legd== TRUE){
                legend(x=1, y=ymax*0.98, legend=colnames(QC), lty=1:ncol(QC), 
                    pch=1:ncol(QC), col=rainbow(ncol(QC)), ncol=2, cex=0.8)
            }
        }     
        title(main=toupper(title), xlab="", ylab="")
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        dev.off()
    }

    ## Start HTML file
    HTMLStart(filename=outputFile)
    titl <- as.title("<p align=center>Quality Assessment</p>")
    HTML(titl, file = outputFile , append = FALSE)
    HTML(targetSort, file = outputFile, Border=1, innerBorder=1)
    HTML("<li> Figure 1 -<a href= #fig1> Intensity Distribution </a>", 
        file = outputFile )
    HTML("<li> Figure 2 -<a href= #fig2> Average Background </a>", 
        file = outputFile )
    HTML("<li> Figure 3 -<a href= #fig3> Scaling Factor </a>", 
        file = outputFile )
    HTML("<li> Figure 4 -<a href= #fig4> Hybridization Controls </a>", 
        file = outputFile )
    HTML("<li> Figure 5 -<a href= #fig5> Housekeeping Controls </a>", 
        file = outputFile )
    HTML("<li> Figure 6 -<a href= #fig6> RNA Quality Assessment </a>", 
        file = outputFile )
    HTML("<li> Figure 7 -<a href= #fig7> Hierarchical Clustering of Samples </a>", 
        file = outputFile )
    HTML("<li> Figure 8 -<a href= #fig8> Pseudo-chip Images </a>", 
        file = outputFile )
    HTML("<li> Figure 9 -<a href= #fig9> Normalized Unscaled Standard Error Plot and Relative Log Expression Plot </a>", 
        file = outputFile )

    ## Figure 1
    figure1 <- "Figure1.png"
    HTML("<hr><a name= \"fig1\"></a>", file = outputFile)
    intensityPlot(figure1, object, targetTag, mycol=rainbow(nchip), mydir=QA.dir)    
    HTMLInsertGraph(figure1, file = outputFile, Caption = "Figure 1: Raw Intensity", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    figure2 <- "Figure2.png"
    figure3 <- "Figure3.png"
    figure4 <- "Figure4.png"
    figure5 <- "Figure5.png"
  
    if (is.na(included)) {
        HTML("<hr><a name= \"fig2\"></a>", file = outputFile)
        emptyPlot(figure2)     
        HTMLInsertGraph(figure2, file = outputFile, 
            Caption = "Figure 2: Average Background/Percentage Present", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

        HTML("<hr><a name= \"fig3\"></a>", file = outputFile)
        emptyPlot(figure3)     
        HTMLInsertGraph(figure3, file = outputFile, 
            Caption = "Figure 3: Scaling Factor", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

        HTML("<hr><a name= \"fig4\"></a>", file = outputFile)
        emptyPlot(figure4)     
        HTMLInsertGraph(figure4, file = outputFile, 
            Caption = "Figure 4: Hybridization Controls", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

        HTML("<hr><a name= \"fig5\"></a>", file = outputFile)
        emptyPlot(figure5)     
        HTMLInsertGraph(figure5, file = outputFile, 
            Caption = "Figure 5: Housekeeping Controls", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    }
    else {
        data_mas5 <- call.exprs(object, sc = 500, "mas5")
        qc <- qc(object, data_mas5)

        ## Figure2
        HTML("<hr><a name= \"fig2\"></a>", file = outputFile)
        png(filename = figure2)
        ##layout(matrix(c(1, 2), nr = 1, ncol = 2), TRUE)
        ##margin=c(10, 4, 3, 1)
        par(mfrow = c(1, 2), las = 2, mar = c(8, 2, 2, 1))
        barplot(avbg(qc), names = targetTag, cex.axis = 0.7, col = rainbow(nchip), 
            cex.main = 0.9, font.main = 2, main = "Average background", 
            las=2)
        barplot(percent.present(qc), names = targetTag, cex.axis = 0.7, 
            col = rainbow(nchip), cex.main = 0.9, font.main = 2, 
            main = "Percentage Present",las=2)
        par(mar=c(6, 4.1, 4.1, 2.1))
        dev.off()
        HTMLInsertGraph(figure2, file = outputFile, 
            Caption = "Figure 2: Average Background/Percentage Present", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=540)
        
        ## Figure3
        HTML("<hr><a name= \"fig3\"></a>", file = outputFile)
        lineGraph(figure3, as.data.frame(sfs(qc)), targetTag, "Scaling Factor", legd=FALSE)
        HTMLInsertGraph(figure3, file = outputFile, Caption = "Figure 3: Scaling Factor", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)
 
        ## Figure4
        HTML("<hr><a name= \"fig4\"></a>", file = outputFile)
        spikeInQC <- data.frame(spikeInProbes(qc))
        colnames(spikeInQC) <- c("BioB", "BioC", "BioD", "Cre")
        lineGraph(figure4, spikeInQC , targetTag, "Hybridization Controls")
        HTMLInsertGraph(figure4, file = outputFile, Caption = "Figure 4: Hybridization Controls", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

        ## Figure5
        HTML("<hr><a name= \"fig5\"></a>", file = outputFile)
        qc2 <- data.frame(2^ratios(qc))
        colnames(qc2) <- c("ACTIN 3'/5'", "GAPDH 3'/5'", "ACTIN 3'/M", "GAPDH 3'/M")
        lineGraph(figure5, qc2 , targetTag, "Housekeeping Controls")
        HTMLInsertGraph(figure5, file = outputFile, Caption = "Figure 5: Housekeeping Controls", 
            GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)
    }

    ## Figure 6
    figure6 <- "Figure6.png"
    deg <- AffyRNAdeg(object)
    flag <- FALSE
    for (i in 1:dim(deg$means.by.number)[1]) {
        for (j in 1:dim(deg$means.by.number)[2]) {
            if (deg$means.by.number[i, j] < 0) {
                deg$means.by.number[i, j] <- NA
                flag <- TRUE
            }
        }
    }
    if (flag == TRUE) {
        for (i in 1:dim(deg$means.by.number)[1]) {
            for (j in 1:dim(deg$means.by.number)[2]) {
                if (is.na(deg$means.by.number[i, j])) {
                  deg$means.by.number[i, j] <- mean(deg$means.by.number[i,], na.rm = TRUE)
                  deg$ses[i, j] <- mean(deg$ses[i, ], na.rm = TRUE)
                }
            }
        }
    }

    png(filename = figure6)
    plotAffyRNAdeg(deg, cols = rainbow(nchip))
    legend("topleft", col = rainbow(nchip), lty = 1, cex = 0.7, legend = targetTag)
    dev.off()
    HTML("<hr><a name= \"fig6\"></a>", file = outputFile)
    HTMLInsertGraph(figure6, file = outputFile, Caption = "Figure 6: RNA Quality Assessment", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=360)

    ## Figure 7
    normal <- rma(object)
    figure7 <- "Figure7.png"
    HTML("<hr><a name= \"fig7\"></a>", file = outputFile)
    plotCluster (figure7, exprs(normal), targetTag , "Hierarchical Clustering of Samples",
        mydir=QA.dir)
    HTMLInsertGraph(figure7, file = outputFile, Caption = "Figure 7: Hierarchical 
        Clustering of Samples", GraphBorder = 1, Align = "center", append = TRUE, 
        WidthHTML=540, HeightHTML=360)

    ## Figure 8
    row <- round(sqrt(nchip))
    for (x in row:nchip) {
        if (((row * x) - nchip) >= 0) {
            col <- x
            break
        }
    }  
    figure8 <- "Figure8.png"
    png(filename = figure8)
    sampleNames(object) <- targetTag 
    data_PLM <- fitPLM(object, output.param = list(varcov = "none"))
    par(mar = c(1, 1, 2, 1), mfrow = c(row, col))
    for (i in 1:nchip) image(data_PLM, which = i)
    dev.off()
    HTML("<hr><a name= \"fig8\"></a>", file = outputFile)
    HTMLInsertGraph(figure8, file = outputFile, Caption = "Figure 8: Pseudo-chip Images", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=540)

    figure9 <- "Figure9.png"
    png(filename = figure9)
    par(mfrow = c(1, 2), las = 2, mar = c(6, 2, 2, 1))
    boxplot(data_PLM, main = "NUSE Plot", cex.axis = 0.7, names = targetTag, col = rainbow(nchip))
    Mbox(data_PLM, main = "RLE Plot", names = targetTag, cex.axis = 0.7, col = rainbow(nchip))
    dev.off()
    HTML("<a name= \"fig9\"></a>", file = outputFile)
    HTMLInsertGraph(figure9, file = outputFile, Caption = "Figure 9: Normalized Unscaled Standard Error Plot and Relative Log Expression Plot", 
        GraphBorder = 1, Align = "center", append = TRUE, WidthHTML=540, HeightHTML=540)
    
    HTML("<hr>This report was generated by ", file = outputFile)
    HTML(paste("<li>", "affy Version", package.version("affy"), 
        "(by Rafael A. Irizarry, Laurent Gautier and Benjamin Bolstad),"), 
        file = outputFile)
    HTML(paste("<li>", "affyPLM Version", package.version("affyPLM"), 
        "(by Benjamin Bolstad),"), file = outputFile)
    HTML(paste("<li>", "limma Version", package.version("limma"), 
        "(by Gordon Smyth),"), file = outputFile)
    HTML(paste("<li>", "R2HTML Version", package.version("R2HTML"), 
        "(by Eric Lecoutre),"), file = outputFile)
    HTML(paste("<li>", "simpleaffy Version", package.version("simpleaffy"), 
        "(by Crispin J. Miller), and"), file = outputFile)
    HTML(paste("<hr>Generated on: ", date()), file = outputFile)
    HTMLStop()

    setwd(cur.dir)

}

