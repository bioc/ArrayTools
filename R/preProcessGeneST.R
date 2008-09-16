`preProcessGeneST` <-
function(object, offset = 1, rmControl = TRUE, output = FALSE, mydir=getwd()){
    
    if(!is.element(annotation(object), c("hugene10st", "mogene10st"))){
        stop("We only support 'hugene10st' or 'mogene10st' for genechip array")
    }   

    if (rmControl) {        
        if (annotation(object) == "hugene10st"){
            data(hugene10stCONTROL)
            control <- hugene10stCONTROL
        } else if (annotation(object) == "mogene10st"){
            data(mogene10stCONTROL)
            control <- mogene10stCONTROL
        }
        exprs <- exprs(object)
        rmIndex <- is.element(rownames(exprs), control$probeset_id)
        exprs <- exprs[!rmIndex,]
        logExprs <- log2(exprs + offset)
    } else {
        logExprs <- log2(exprs(object) + offset)
    }
    exprs(object) <- logExprs
    if (output) {
        dir <- mydir
        packageName <- paste(annotation(object), ".db", sep="")
        library(packageName , character.only=TRUE)
        require("annaffy")      
        probeset_id = rownames(logExprs)
        Symbol <- getText(aafSymbol(probeset_id, packageName))
        Description <- getText(aafDescription(probeset_id, packageName))
        Chromosome <- getText(aafChromosome(probeset_id, packageName))
        GenBank <- getText(aafGenBank(probeset_id, packageName))
        Cytoband.temp <- lapply(aafCytoband(probeset_id, packageName), function(x) x@band)
        Cytoband.temp1 <- lapply(Cytoband.temp, function(x){
           if (length(x) == 0) return ("") else return (x)})
        Cytoband <- sapply(Cytoband.temp1, paste, collapse = "; ")
        UniGene <- getText(aafUniGene(probeset_id, packageName))
        PubMed <- getText(aafPubMed(probeset_id, packageName))
        LocusLink <- getText(aafLocusLink(probeset_id, packageName))   
        logExprs <- data.frame(probeset_id = probeset_id, logExprs,
           Symbol = Symbol, Description = Description, Chromosome = Chromosome,
           GenBank = GenBank, Cytoband = Cytoband, UniGene = UniGene, 
           PubMed = PubMed, LocusLink = LocusLink )
               
        result.dir <- paste(dir, "Result", sep="/")
        dir.create(path=result.dir, showWarnings = FALSE)
        setwd(result.dir)
        normal.dir <- paste(result.dir, "Normalized_Data", sep="/")
        dir.create(path=normal.dir, showWarnings = FALSE)
        setwd(normal.dir)
        write.csv(logExprs, file="normaldata.csv", row.names=F)
        setwd(dir)
    }
    experimentData(object)@preprocessing <- list(offset = offset, 
        rmControl = rmControl)
    return(object)
}

