
`output.ing` <-
function (allfile, eSet, filename = "IngenuityFile") 
{
    data1 <- allfile
    dat1 <- data1[[1]]
    filename1 <- getFileName(dat1)
    if (length(filename1) > 1) {
        filename1 <- filename1[[1]]
        result1 <- regress(object = eSet, contrast = getContrast(dat1)[[1]], 
            method = regressionMethod(dat1)[[1]], adj = adjustment(dat1))[[1]]
    }
    else {
        result1 <- regress(object = eSet, contrast = getContrast(dat1), 
            method = regressionMethod(dat1), adj = adjustment(dat1))
    }
    FC1 <- as.data.frame(getFC(result1))
    if (ncol(FC1) > 1) 
        stop("We can only create an ingenuity file for one comparison")
    top.Table1 <- data.frame(getID(result1), FC1, getF(result1), 
        getP(result1), getAdjP(result1))
    colnames(top.Table1) <- c("probeset_id", paste(filename1, 
        c("logFC", "F", "P", "AdjP"), sep = ":"))
    if (length(data1) > 1) {
        for (i in 1:(length(data1) - 1)) {
            dat <- data1[[i + 1]]
            filename <- getFileName(dat)
            if (length(filename) > 1) {
                filename <- filename[[1]]
                result <- regress(object = eSet, contrast = getContrast(dat)[[1]], 
                  method = regressionMethod(dat)[[1]], adj = adjustment(dat)[[1]])
            }
            else {
                result <- regress(object = eSet, contrast = getContrast(dat), 
                  method = regressionMethod(dat), adj = adjustment(dat))
            }
            FC <- as.data.frame(getFC(result))
            if (ncol(FC) > 1) 
                stop("We can only create an ingenuity file for one \n                
                    comparison")
            top.Table <- data.frame(FC, getF(result), getP(result), 
                getAdjP(result))
            colnames(top.Table) <- paste(filename, c("logFC", 
                "F", "P", "AdjP"), sep = ":")
            top.Table1 <- cbind(top.Table1, top.Table)
        }
    }
    packageName <- paste(getAnnotation(allfile[[1]]), ".db", sep="")
    library(packageName , character.only=TRUE)
    require("annaffy")  
    probeset_id = as.character(top.Table1$probeset_id)
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
    file.ann <- data.frame(top.Table1, Symbol = Symbol, Description = Description, 
        Chromosome = Chromosome, GenBank = GenBank, Cytoband = Cytoband, 
        UniGene = UniGene, PubMed = PubMed, LocusLink = LocusLink )
    write.table(file.ann, "Ingenuity.txt", sep = "\t", row.names = F, 
        quote = F)
}


