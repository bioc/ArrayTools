
`createIndex` <- 
function(..., mydir=getwd(), index.file= "index.html", 
    createHeader = NULL){    

    result <- list(...)  
   
    ##browser()
    library(xtable)
    arrayType <- getAnnotation(result[[1]])
    if (!all(sapply(result, getAnnotation) == arrayType))
        stop ("Array types have to be the same for each object!")

    if (is.element(arrayType, c("hugene10st", "mogene10st"))) {  ## For genechip array only
        offsetList <- sapply(lapply(result, getNormalizationMethod), function(x) x$offset)
        offset <- offsetList[1]
        rmControlList <- sapply(lapply(result, getNormalizationMethod), function(x) x$rmControl)
        rmControl <- rmControlList[1]
        if (!all(offsetList == offset) && all(rmControlList == rmControl))
            stop ("Normalization methods have to be the same for each object!")
    }  else {
        normalMethodList <- sapply(lapply(result, getNormalizationMethod), function(x) x$normalMethod)
        normalMethod <- normalMethodList[1]
        if (!all(normalMethodList == normalMethod))
            stop ("Normalization methods have to be the same for each object!")
    }
    
    ## Get Filtering Information
    numChipList <- sapply(lapply(result, getFilterMethod), function(x) x$numChip)
    numChip <- numChipList[1]
    bgList <- sapply(lapply(result, getFilterMethod), function(x) x$bg)
    bg <- bgList[1]
    rangeList <- sapply(lapply(result, getFilterMethod), function(x) x$range)
    range <- rangeList[1]
    iqrPctList <- sapply(lapply(result, getFilterMethod), function(x) x$iqrPct)
    iqrPct <- iqrPctList[1]
    if (!all(numChipList == numChip) && all(bgList == bg) &&
        all(rangeList == range) && all(iqrPctList == iqrPct))
        stop ("Filtering methods have to be the same for each object!")

    numResult <- length(result)
    regressResultIndex <- (1:numResult)[is.element(sapply(result, class), "regressResult")]
    interactionResultIndex <- (1:numResult)[is.element(sapply(result, class), "interactionResult")]
    
    if (length(regressResultIndex)>0) 
        regressResult <- result[regressResultIndex]
    else 
        regressResult <- NULL

    if (length(interactionResultIndex)>0) 
        interactionResult <- result[interactionResultIndex]
    else 
        interactionResult <- NULL
  
    curDir <- mydir
    setwd(curDir)
    result.dir <- paste(curDir, "Result", sep="/")
    if (!file.exists("Result")) stop ("Your current directory doesn't contain the 'Result' folder")
    setwd(result.dir)  

    ## Start HTML file
    HTMLStart(filename=index.file)
    titl <- as.title("<p align=left><font size='6', face='Verdana'>Analysis Report</p>")
    HTML(titl, file = index.file, append=FALSE)

    if (!is.null(createHeader)){          
        HTML(paste("<p align=left><font size='2', face='Oxford'>", 
            paste(createHeader, collapse=" <BR> "), sep= " "), file = index.file)
    }
    
    if (length(getContrast(result[[1]])) > 1){
        target <- getTarget(getContrast(result[[1]][[1]]))
    } else {
        target <- getTarget(getContrast(result[[1]]))
    }

    ####################Insert target file####################
    HTML("<hr /> <p align=left><font size='5', face='Verdana'> Sample", file = index.file)
    HTML("<ul>", file=index.file)
    HTML(target, file=index.file, innerBorder=2, align="left")
    HTML("</ul>", file=index.file)
   

    ##################Insert 'Data Preprocessing' information##################
    if (file.exists("Normalized_Data")) {    
        HTML("<hr /> <p align=left><font size='5', face='Arial'> Data Preprocessing", 
            file = index.file)
        HTML("<ul>", file=index.file)     
        if (is.element(arrayType, c("hugene10st", "mogene10st"))){            
            HTML("<p align=left><li><font size='3', face='Arial'> All data processing was done 
                using Bioconductor or APT tools", file = index.file)
            if (offset > 0){
                HTML(paste("<p align=left><li><font size='3', face='Arial'> The offset of", 
				offset, "was added to the origianl value before taking the log2", sep=" "), 
				file = index.file) 
            }
            if (rmControl){
                HTML("<p align=left><li><font size='3', face='Arial'> 'normgene' and 'control'
                    genes were removed", file = index.file)
            }
        } else {
            HTML(paste("<p align=left><li><font size='3', face='Arial'> Background Correction, normalization,
                and summarization used", toupper(normalMethod), "method", sep=" "), file = index.file)
        }
        if (!is.null(regressResult)){
            nSample <- getFilterMethod(regressResult[[1]])$nSample
            nGene <- getFilterMethod(regressResult[[1]])$nGene
        } else {
            nSample <- getFilterMethod(interactionResult[[1]][[1]])$nSample
            nGene <- getFilterMethod(interactionResult[[1]][[1]])$nGene
        }

        HTML(paste("<p align=left><li><font size='3', face='Arial'> Endpoint: data matrix of", 
            nGene, "genes X", nSample, "samples", sep=" "),file = index.file)
        HTML("<p align=left><li><font size='3', face='Arial'> <a href=Normalized_Data\\normaldata.csv> 
           Download Processed Data </a>", file = index.file)
        HTML("</ul>", file=index.file)

    }

    ###########################Insert 'QA' link###########################

    if (file.exists("Quality_Control_Assessment")) {
 
        HTML("<hr /> <p align=left><font size='5', face='Arial'> Quality Assessment", file = index.file)
        HTML("<ul>", file=index.file)

        if (is.element(arrayType, c("hugene10st", "mogene10st"))){
            HTML("<p align=left> <font size='3'> <a href=Quality_Control_Assessment\\QA.html> Quality Assessment 
                Report </a> includes the following figures: Intensity Distribution, Mean Signal, BAC Spike, 
                polya Spike, Pos Vs Neg Auc, Mad Residual Signal, RLE MEAN, and Hierarchical Clustering of Samples
                plots", file = index.file)
            
        } else{
            HTML("<p align=left> <font size='3'> <a href=Quality_Control_Assessment\\QA.html> Quality Assessment 
                Report </a> includes the following figures: Intensity Distribution, Average Background, Scaling Factor, 
                Hybridization Controls, Housekeeping Controls, RNA Quality Assessment, 
                Hierarchical Clustering of Samples, Pseudo-chip Images, and 
                Normalized Unscaled Standard Error Plot and Relative Log Expression Plots", file = index.file)
        }
        HTML("</ul>", file=index.file)

	}

    #######################Insert 'Filtering' link#######################

    if (file.exists("Filtered_Data")) {

        HTML("<hr /> <p align=left><font size='5', face='Arial'> Filtering", file = index.file)
        HTML("<ul>", file=index.file)

        if (numChip > 0 && bg > 0) {
        
            HTML(paste("<p align=left><li><font size='3', face='Arial'> Probes with log2 intensity <", 
                bg, "in more than", nSample - numChip, "arrays are considered uninformative 
                and removed", sep=" "), file = index.file)
        }
        
        if (range > 0 ) {
            HTML(paste("<p align=left><li><font size='3', face='Arial'> Probes with range <", 
                range, "are considered uninformative and removed", sep=" "), file = index.file)
        }

        if (iqrPct > 0 ) {
            Pct <- paste(iqrPct*100, "%")
            HTML(paste("<p align=left><li><font size='3', face='Arial'> Probes with IQR <", 
                Pct, "are considered uninformative and removed", sep=" "), file = index.file)
        }

        if (!is.null(regressResult)){
            nGeneFilter <- length(regressResult[[1]]@ID)
        } else {
            nGeneFilter <- length(interactionResult[[1]][[1]]@ID)
        }

        HTML(paste("<p align=left><li><font size='3', face='Arial'> After filtering:", nGeneFilter, 
           "probesets are left", sep=" "), file = index.file)
        HTML("<p align=left><li><font size='3', face='Arial'> <a href=Filtered_Data\\filterData.csv> 
           Download Filtered Data </a>", file = index.file)
        HTML("</ul>", file=index.file)
    
    }

    ###################Insert 'Diffrerential Expressed Genes'link#########
    if (file.exists("Differentially_Expressed_Genes")){
       
        HTML("<hr /> <p align=left><font size='5', face='Arial'> Identification of Differentially 
            Expressed Genes", file = index.file)        

        if (!is.null(regressResult)) {

            HTML("<ul>", file=index.file)
            covariateList <- lapply(regressResult, function(x) getCovariates(getContrast(x)))
            mainEffect <- sapply(covariateList , function(x) x[1])
            Levels <- sapply(regressResult, function(x) getLevel(getContrast(x)))
            Levels[is.na(Levels)] <- "none"
            IntIndex <- lapply(regressResult, function(x) getIntIndex(getContrast(x)))
            
            adjusting1 <- c()
            for (i in 1:length(IntIndex)){
                if (identical(IntIndex[[i]], 0)) adjusting1[[i]] <- covariateList [[i]][-1]
                else adjusting1[[i]] <- covariateList [[i]][-IntIndex[[i]]]
            }

            Adjusting <- sapply(adjusting1, paste, collapse=",")
		    Adjusting[Adjusting==""] <- "none"
            
            compare1 <- sapply(regressResult, function(x) getCompare1(getContrast(x)))
            compare2 <- sapply(regressResult, function(x) getCompare2(getContrast(x)))
            Comparison <- paste(compare1, compare2, sep=" vs. ")
            
            Ratio <- paste("log2(", 2^(sapply(regressResult, function(x) getFCCutoff(x))), ")", sep="")
            P  <- sapply(regressResult, function(x) getPCutoff(x))
            MultCompAdjtmt <- sapply(regressResult, function(x) adjustment(x))
            StatMethod <- sapply(regressResult, function(x) regressionMethod(x))
        
            link <- character(length(regressResult))

            for (i in 1:length(regressResult)){
                fname <- paste(getFileName(regressResult[[i]]), "html", sep=".")
                nSigGene <- sum(getIndex(regressResult[[i]]))
                link[i] <- paste("<a href=Differentially_Expressed_Genes\\", fname, "> ", nSigGene, 
                    " genes", "</a>", sep="")
            }

            withAdjust <- !(sum(is.element(Adjusting, "none")) == length(Adjusting))
            withLevel <- !(sum(is.element(Levels, "none")) == length(Levels))
            
            if (withAdjust && withLevel) 
                resultTable1 <- data.frame(mainEffect=mainEffect, Comparison =Comparison,
                    Levels = Levels, Adjusting = Adjusting, Ratio = Ratio, P = P, 
                    MultCompAdjtmt = MultCompAdjtmt, StatMethod = StatMethod,  Results=link)
            else if (withAdjust && !withLevel)
                resultTable1 <- data.frame(mainEffect=mainEffect, Comparison =Comparison,
                    Adjusting = Adjusting, Ratio = Ratio, P = P, 
                    MultCompAdjtmt = MultCompAdjtmt, StatMethod = StatMethod,  Results=link)
            else if (!withAdjust && withLevel) 
                resultTable1 <- data.frame(mainEffect=mainEffect, Comparison =Comparison,
                    Levels = Levels, Ratio = Ratio, P = P, 
                    MultCompAdjtmt = MultCompAdjtmt, StatMethod = StatMethod,  Results=link)
            else 
	          resultTable1 <- data.frame(mainEffect=mainEffect, Comparison =Comparison,
                Ratio = Ratio, P = P, MultCompAdjtmt = MultCompAdjtmt, 
                StatMethod = StatMethod,  Results=link)
            
            resultTable2 <- xtable(resultTable1)
            print(resultTable2, type="html", file=index.file, include.rownames = FALSE, 
                sanitize.text.function = function(x) x, append=T)
            HTML("</ul>", file=index.file)
        } 

        if (!is.null(interactionResult)) {

            HTML("<ul>", file=index.file)
            for (i in 1:length(interactionResult)){
                covariateList <- getCovariates(getContrast(interactionResult[[i]])[[1]])
                Levels <- c("none", sapply(getContrast(interactionResult[[i]]), getLevel)[-c(1,2)])
		    IntIndex <- getIntIndex(getContrast(interactionResult[[i]])[[1]])
                adjusting1 <- rep(covariateList[-IntIndex], length(Levels))
                Adjusting <- sapply(adjusting1, paste, collapse=",")

                compare1 <- sapply(getContrast(interactionResult[[i]]), getCompare1)[-c(1,2)]
                compare2 <- sapply(getContrast(interactionResult[[i]]), getCompare2)[-c(1,2)]
		        Comparison <- c(paste(covariateList[IntIndex], collapse = " X "),
                    paste(compare1, compare2, sep=" vs. "))
                Ratio <- paste("log2(", 2^unlist(getFCCutoff(interactionResult[[i]]))[-2], ")", sep="")
                Ratio[Ratio == "log2(1)"] <- "none"
                P <- unlist(getPCutoff(interactionResult[[i]]))[-2]
	            MultCompAdjtmt <- unlist(adjustment(interactionResult[[i]]))[-2]
                StatMethod  <- unlist(regressionMethod(interactionResult[[i]]))[-2]
                Results <- sapply(getIndex(interactionResult[[i]]), sum)[-2]
              
                fName <- getFileName(interactionResult[[i]])[[1]]
                
                HTML(paste("<p align=left><font size='3'> <a href=Differentially_Expressed_Genes\\",
                     fName, ".html>", fName, "</a>", sep=""), file = index.file)   
                  
                if (length(Adjusting) > 0) 
                    resultTable1 <- data.frame(Comparison =Comparison,
                        Levels = Levels, Adjusting = Adjusting, Ratio = Ratio, P = P, 
                        MultCompAdjtmt = MultCompAdjtmt, StatMethod = StatMethod,  Results=Results)
                else 
			    resultTable1 <- data.frame(Comparison =Comparison,
                        Levels = Levels, Ratio = Ratio, P = P, 
                        MultCompAdjtmt = MultCompAdjtmt, StatMethod = StatMethod,  Results=Results)

	          resultTable2 <- xtable(resultTable1)
                print(resultTable2, type="html", file=index.file, include.rownames = FALSE, 
                     sanitize.text.function = function(x) x, append=T)
                
            }        
            HTML("</ul>", file=index.file)   
        }
    }

    if (file.exists("Pathway_Analysis")){

        #############################Insert 'Pathway'link#############################
        HTML("<hr /> <p align=left><font size='5', face='Arial'> Pathway Analysis",  file = index.file)
        HTML(paste("<li><p align=left><font size='4', face='Arial'><a href= Pathway_Analysis\\GSEA >",
            "Gene Set Enrichment Analysis",  "</a>", sep=""), file = index.file)
    
        HTML("</ul>", file=index.file)
        HTML(paste("<li><p align=left><font size='4', face='Arial'><a href= http://www.ingenuity.com >", 
            "Ingenuity Pathway Analysis", "</a>", sep=""), file = index.file)
    }

    HTMLStop()
    setwd(curDir)

}
