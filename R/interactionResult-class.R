###############  interactionResult Class ##################
setClass("interactionResult", contains="list",
   validity=function(object) {
       ok <- sapply(slot(object, ".Data"), is, "regressResult")
       if (!all(ok)) "values must be regressResult"
       else TRUE
})
setMethod("getLength", signature("interactionResult"), function(object) length(object@.Data))
setMethod("getID", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@ID))
setMethod("getFC", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@foldChange))
setMethod("getF", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@FValue))
setMethod("getP", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@pValue))
setMethod("getAdjP", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@adjPVal))
setMethod("getContrast", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@contrast))
setMethod("regressionMethod", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@regressionMethod))
setMethod("adjustment", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@adjustment))
setMethod("getIndex", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@significantIndex))
setMethod("getPCutoff", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@significantPvalueCutoff))
setMethod("getFCCutoff", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@significantFCCutoff))
setMethod("getFileName", signature("interactionResult"), function(object)
    lapply(object, function (object) object@fileName))
setMethod("getAnnotation", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@annotation)[[1]])
setMethod("getNormalizationMethod", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@normalizationMethod)[[1]])
setMethod("getFilterMethod", signature("interactionResult"), function(object) 
    lapply(object, function (object) object@filterMethod)[[1]])

setMethod("Sort", signature("interactionResult"), 
    function(x, sorted.by = c("pValue", "log2Ratio", "F"), top=20){

    sorted.by<- match.arg(sorted.by) 

    cat("All genes: \n")
    allTable <- Sort(x[[1]], sorted.by = sorted.by, top = top)
    cat("\nGenes without interaction:\n")
    TableWithoutInt <- Sort(x[[2]], sorted.by = sorted.by, top = top)
    tab1 <- c()
    for (i in 1:(length(x)-2)){
        cat("\nGenes with interaction, level = ", getLevel(getContrast(x[[i+2]])), "\n")
        tab1 <- c(tab1, list(Sort(x[[i+2]], sorted.by = sorted.by, top = top)))
    }
    tab <- c(list(allTable), list(TableWithoutInt), tab1)
    invisible(tab)
})
setMethod("summary", signature("interactionResult"), 
    function(object, ...){

    ##browser()
    cat("All genes: N =", length(getID(object[[1]])), "\n")
    allTable <- summary(object[[1]], ...)
    cat("\nGenes without interaction: N =", length(getID(object[[2]])), "\n")
    TableWithoutInt <- summary(object[[2]], ...)
    tab1 <- c()
    for (i in 1:(length(object)-2)){
        cat("\nGenes with interaction, level = ", getLevel(getContrast(object[[i+2]])), 
        ": N = ", length(getID(object[[i+2]])), "\n")
        tab1 <- c(tab1, list(summary(object[[i+2]], ...)))
    }
    tab <- c(list(allTable), list(TableWithoutInt), tab1)
    invisible(tab)
})
setMethod("Output2HTML", signature("interactionResult"), 
    function (object, mydir=getwd())
{   
    intIndex <- apply(as.data.frame(getIndex(object)[-c(1,2)]), 1, sum) >0    
    if(sum(intIndex) == 0) stop("No significant result is generated!")


    dir <- mydir
    result.dir <- paste(dir, "Result", sep="/")
    dir.create(path=result.dir,  showWarnings = FALSE)
    setwd(result.dir)
    
    diff.dir <- paste(result.dir, "Differentially_Expressed_Genes", sep="/")
    dir.create(path=diff.dir, showWarnings = FALSE)
    setwd(diff.dir)

    packageName <- paste(getAnnotation(object), ".db", sep="")
    library(packageName , character.only=TRUE)
    require("annaffy")    
 
    probeset_id <- getID(object)[[3]][intIndex]
    FCStrata <- as.data.frame(getFC(object)[-c(1,2)])
    FCStrata.select <- FCStrata[intIndex,, drop=F]

    PStrata <- as.data.frame(getP(object)[-c(1,2)])
    PStrata.select <- PStrata[intIndex,, drop=F]

    AdjPStrata <- as.data.frame(getAdjP(object)[-c(1,2)])
    AdjPStrata.select <- AdjPStrata[intIndex,, drop=F]

    AllData <- c()
    for(i in 1:ncol(FCStrata)){
        AllData [[3*i-2]] <- FCStrata.select[[i]]
        AllData [[3*i-1]] <- PStrata.select[[i]] 
        AllData [[3*i]] <- AdjPStrata.select[[i]]
    }  

    Levels <- sapply(getContrast(object), getLevel)[-c(1,2)]
    FCNames <- paste(Levels, "FC", sep=":")
    PNames <- paste(Levels, "P", sep=":")
    adjPNames <- paste(Levels, "adjP", sep=":")
    AllNames <-c()
    for(i in 1:ncol(FCStrata)){
        AllNames[3*i-2] <- FCNames[i]
        AllNames[3*i-1] <- PNames[i] 
        AllNames[3*i] <- adjPNames[i]
    }
    names(AllData) <- AllNames
    AllData <- as.data.frame(AllData)

    IDIndex <- is.element(getID(object[[1]]), probeset_id)
    PInt <- as.data.frame(getP(object)[[1]])
    PInt.select <- PInt[IDIndex,, drop=F]
    names(PInt.select) <- "Int.P"
    AdjPInt <- as.data.frame(getAdjP(object)[[1]])
    AdjPInt.select <- AdjPInt[IDIndex,, drop=F]
    names(AdjPInt.select) <- "Int.AdjP"

    if (sum(is.na(PInt)) == length(getP(object)[[1]])){
        fullresult <- data.frame(Probeset_id  = probeset_id, AllData[c(1,4)])
    } else {
        fullresult <- data.frame(Probeset_id  = probeset_id, AllData,
            PInt.select, AdjPInt.select)       
    }
    anncols <- aaf.handler()[c(2,3,4,6,7,8,9,10)]
    anntable <- aafTableAnn(probeset_id, packageName , anncols)
    testtable <- aafTable(items=fullresult)
    table <- merge(testtable, anntable)
    saveHTML(table, paste(getFileName(object[[1]]), "html", sep = "."), title = "Target Genes")
    setwd(dir)
})
setMethod("show", signature("interactionResult"), 
    function(object){ 

    cat("All genes: \n")
    allTable <- show(object[[1]])
    cat("\nGenes without interaction:\n")
    TableWithoutInt <- show(object[[2]])
    tab1 <- c()
    for (i in 1:(length(object)-2)){
        cat("\nGenes with interaction, level = ", getLevel(getContrast(object[[i+2]])), "\n")
        tab1 <- c(tab1, list(show(object[[i+2]])))
    }
    tab <- c(list(allTable), list(TableWithoutInt), tab1)
    invisible(tab)

})


