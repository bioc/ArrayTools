`selectSigGene` <-
function(object, p.value = 0.05, fc.value = 0) {    
    
    if (sum(is.na(getP(object))) == length(getP(object)))
        good.p <- is.na(getP(object))
    else {
  if (adjustment(object) == "none")
            good.p <- getP(object) < p.value
        else
            good.p <-getAdjP(object) < p.value
    }

    FC <- as.data.frame(getFC(object))
    good.matrix <- matrix(0, nrow = nrow(FC), ncol = ncol(FC))
    for (i in 1:ncol(FC)) {
        good.matrix[, i] <- abs(FC[,i]) >= fc.value
    }
    good.fc <- (apply(good.matrix, 1, sum) > 0)

    object@significantIndex <- good.p & good.fc
    object@significantPvalueCutoff <- p.value
    object@significantFCCutoff <- fc.value

    varList <- getCovariates(getContrast(object))
    compare1 <- getCompare1(getContrast(object))
    compare2 <- getCompare2(getContrast(object))
    
    if (length(getIntIndex(getContrast(object))) > 1) {
        if(!is.na(compare1)){
            index <- logical(length(varList))
            for (i in 1:length(varList)){
                index[i] <- is.element(compare1, 
                    as.matrix(getTarget(getContrast(object))[varList[i]]))
            }
            mainEffect <- varList[index]
            intVar1 <- intVar2 <-NA
            
        } else{            
            mainEffect <- NA
            intVar1 <- varList[getIntIndex(getContrast(object))[1]]
            intVar2 <- varList[getIntIndex(getContrast(object))[2]]    
        }        
    } else {
        mainEffect <- varList[1]
        intVar1 <- NA
        intVar2 <- NA
    }  
    if (length(varList) > 1)
        if (getInteraction(getContrast(object))) {
            adjVar <- paste(varList[-getIntIndex(getContrast(object))], 
                collapse="_")
            if (adjVar == "") adjVar <- NA
  } else {
            if(is.na(getLevel(getContrast(object))))
                adjVar <- paste(varList[-1], collapse="_")
            else 
                adjVar <- paste(varList[-getIntIndex(getContrast(object))], 
                    collapse="_")
            if (adjVar == "") adjVar <- NA
        }
    else 
        adjVar <- NA

    if(is.na(adjVar)){
        if (getInteraction(getContrast(object))) {
           object@fileName <- paste(intVar1, "-X-", intVar2, sep="")
        }  else {
            if (is.na(getLevel(getContrast(object))))
                object@fileName <- paste(mainEffect, "-", compare1, 
                    ".VS.", compare2, sep="")
            else
     object@fileName <- paste(mainEffect, "-", compare1, 
                    ".VS.", compare2, "-At-", getLevel(getContrast(object)), sep="")
        }
    } else {
        if (getInteraction(getContrast(object))) {
            object@fileName <- paste(intVar1, "-X-", intVar2, "-",
                "ADJ", "-", adjVar, sep="")
        } else {
            if (is.na(getLevel(getContrast(object))))
                object@fileName <- paste(mainEffect, "-", compare1, ".VS.", compare2, 
                    "-", "ADJ", "-", adjVar, sep="")
            else
    object@fileName <- paste(mainEffect, "-", compare1, ".VS.", compare2, 
                    "-", "ADJ", "-", adjVar, "-At-", 
                    getLevel(getContrast(object)), sep="")
        }
    }
    return(object)
}

