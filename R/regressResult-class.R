
###############  RegressResult Class ##################
setClass("regressResult", representation (ID = "character",
    foldChange = "list", FValue = "numeric", pValue = "numeric", 
    adjPVal = "numeric", contrast ="contrastMatrix",
    regressionMethod = "character", adjustment = "character", 
    significantIndex = "logical", significantPvalueCutoff = "numeric",
    significantFCCutoff = "numeric", fileName ="character", annotation="character",
	normalizationMethod = "list", filterMethod = "list"))

setMethod("getID", signature("regressResult"), function(object) object@ID)
setMethod("getFC", signature("regressResult"), function(object) object@foldChange)
setMethod("getF", signature("regressResult"), function(object) object@FValue)
setMethod("getP", signature("regressResult"), function(object) object@pValue)
setMethod("getAdjP", signature("regressResult"), function(object) object@adjPVal)
setMethod("getContrast", signature("regressResult"), function(object) object@contrast)
setMethod("regressionMethod", signature("regressResult"), function(object) object@regressionMethod)
setMethod("adjustment", signature("regressResult"), function(object) object@adjustment)
setMethod("getIndex", signature("regressResult"), function(object) object@significantIndex)
setMethod("getPCutoff", signature("regressResult"), function(object) object@significantPvalueCutoff)
setMethod("getFCCutoff", signature("regressResult"), function(object) object@significantFCCutoff)
setMethod("getFileName", signature("regressResult"), function(object) object@fileName)
setMethod("getAnnotation", signature("regressResult"), function(object) object@annotation)
setMethod("getNormalizationMethod", signature("regressResult"), function(object) object@normalizationMethod)
setMethod("getFilterMethod", signature("regressResult"), function(object) object@filterMethod)

setMethod("Sort", signature("regressResult"), 
    function(x, sorted.by = c("pValue", "log2Ratio", "F"), top=20){

    sorted.by<- match.arg(sorted.by) 
    FC <- as.data.frame(getFC(x))    
    dataFrame <- data.frame(ID=getID(x), FC, F=getF(x), pValue=getP(x), 
        adjPVal=getAdjP(x))
	if (length(getIndex(x)>0)){
        cat("\nThere are ", sum(getIndex(x)), " significant genes.\n\n")
		dataFrame <- dataFrame[getIndex(x),]
    }	
    if (sorted.by == "pValue") {
        if (adjustment(x) == "none") {
            top.Table <- dataFrame[order(dataFrame$pValue),]
        } else {
            top.Table <- dataFrame[order(dataFrame$adjPVal),]
        }
    } else if (sorted.by == "F"){
        top.Table <- dataFrame[order(dataFrame$F, decreasing = T),]
    } else if (sorted.by == "log2Ratio"){
        top.Table <- dataFrame[order(abs(dataFrame[,2]), decreasing = T),]
    }
	if (top > nrow(dataFrame)) top <- nrow(dataFrame)
	Top <- top.Table[1:top, ]
    print(Top)
    
    invisible(Top)
})
setMethod("show", signature("regressResult"), 
    function(object){ 
    FC <- as.data.frame(getFC(object)) 
    if (length(getIndex (object)) > 0)  
        dataFrame <- data.frame(ID=getID(object), FC, F=getF(object), 
            pValue=getP(object), adjPVal=getAdjP(object), significant=getIndex (object))
    else 
	    dataFrame <- data.frame(ID=getID(object), FC, F=getF(object), 
            pValue=getP(object), adjPVal=getAdjP(object))        

    print(dataFrame[1:20,])
    invisible(dataFrame)
})
setMethod("summary", signature("regressResult"), 
    function(object, option=c("cummulative", "discrete"), ...){
    ##browser()
    option<- match.arg(option)
   
    cumTable <- function(tab) {
        ##browser()
        ##index <- which(rownames(tab) == 'H:(-log2(1.3), 0]')   

        temp <- is.element(substr(rownames(tab), 1, 1), LETTERS[1:8])
        index <- max(which(temp ==TRUE))
     
        tab1 <- apply(as.matrix(tab[index:1,]), 2, cumsum)
        tab2 <- apply(as.matrix(tab[(index+1):nrow(tab),]), 2, cumsum)
        if (ncol(tab1) ==1) {
            tab1 <- t(tab1)
            rownames(tab1) <- rownames(tab)[index:1]
        }
        if (ncol(tab2) ==1) {
            tab2 <- t(tab2)
            rownames(tab2) <- rownames(tab)[(index+1):nrow(tab)]
        }
        newtab <- rbind(tab1[index:1,, drop=F], tab2)
    }
    P <- cut(getP(object), breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), right = F, 
        labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">=0.05"))
    FC <- as.data.frame(getFC(object)) 
    cut1 <- c(-Inf, -log2(c(5, 4, 3, 2.5, 2, 1.5, 1.3)), 
        log2(c(1, 1.3, 1.5, 2, 2.5, 3, 4, 5)), Inf)

    lab <- c("A:(-Inf, -log2(5)]", "B:(-log2(5), -log2(4)]", "C:(-log2(4), -log2(3)]",
        "D:(-log2(3), -log2(2.5)]", "E:(-log2(2.5), -log2(2)]", "F:(-log2(2), -log2(1.5)]",
        "G:(-log2(1.5), -log2(1.3)]", "H:(-log2(1.3), 0]", "I:(0, log2(1.3)]",
        "J:(log2(1.3), log2(1.5)]", "K:(log2(1.5), log2(2)]", "L:(log2(2), log2(2.5)]",
        "M:(log2(2.5), log2(3)]", "N:(log2(3), log2(4)]", "O:(log2(4), log2(5)]", "P:(log2(5), Inf)")

    FoldChange <- apply(FC, 2, function(x) return(cut(x, breaks = cut1, include.lowest = T, labels=lab)))
    nFC <- ncol(FoldChange)
    if (nFC > 1) {
        tab <- c()
        for (i in 1:nFC){
            tabTemp <- table(FoldChange[,i], P) 
            if (option == "cummulative") tabTemp <- cumTable(tabTemp)  
            rownames(tabTemp) <- substr(rownames(tabTemp), 3, sapply(rownames(tabTemp), nchar))
            tab <- c(tab, list(tabTemp))
        }
        names(tab) <- names(FC)
    } else {
        tab <- table(FoldChange, P)
        if (option == "cummulative") tab <- cumTable(tab)
        rownames(tab) <- substr(rownames(tab), 3, sapply(rownames(tab), nchar))
    }   
    print(tab)
    invisible(tab)
 
})
setMethod("Output2HTML", signature("regressResult"), 
    function (object, mydir=getwd()){    
    
    ##browser()
    if( sum(getIndex(object)) == 0) stop("No significant result is generated!")

    dir <- mydir
    result.dir <- paste(dir, "Result", sep="/")
    dir.create(path=result.dir, showWarnings = FALSE)
    setwd(result.dir)
    
    diff.dir <- paste(result.dir, "Differentially_Expressed_Genes", sep="/")
    dir.create(path=diff.dir, showWarnings = FALSE)
    setwd(diff.dir)

    packageName <- paste(getAnnotation(object), ".db", sep="")
    library(packageName , character.only=TRUE)
    require("annaffy")      
     
    probeset_id  <- getID(object)[getIndex(object)]
    FC <- as.data.frame(getFC(object))
    FC.select <- FC[getIndex(object),, drop=F]
    P.Value <- getP(object)[getIndex(object)]
    Adj.P.Value <- getAdjP(object)[getIndex(object)]

    if (sum(is.na(P.Value)) == length(P.Value)){
        fullresult <- data.frame(Probeset_id  = probeset_id , FC.select)
    } else {
        fullresult <- data.frame(Probeset_id  = probeset_id , FC.select,
            P.Value =  P.Value, Adj.P.Value = Adj.P.Value)
    }
    ##gN <- as.character(result$ID)
    anncols <- aaf.handler()[c(2,3,4,6,7,8,9,10)]
    anntable <- aafTableAnn(probeset_id, packageName , anncols)
    testtable <- aafTable(items=fullresult)
    table <- merge(testtable, anntable)
    saveHTML(table, paste(getFileName(object), "html", sep = "."), title = "Target Genes")
    setwd(dir)
})


