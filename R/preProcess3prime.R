`preProcess3prime` <-
function (object, method=c("rma", "gcrma"), output = FALSE, mydir=getwd()) 
{
    method <- match.arg(method)
    if (method == "gcrma") {
        require("gcrma")
        normaldata <- gcrma(object)
    } else if (method == "rma") {
        normaldata <- rma(object)
    } else stop("Only 'rma' or 'gcrma' are supported")
    if (output) {
        dir <- mydir
        packageName <- paste(annotation(object), ".db", sep="")
        library(packageName , character.only=TRUE)
        require("annaffy")      
        probeset_id = rownames(exprs(normaldata ))
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
        Exprs <- data.frame(probeset_id = probeset_id, exprs(normaldata),
           Symbol = Symbol, Description = Description, Chromosome = Chromosome,
           GenBank = GenBank, Cytoband = Cytoband, UniGene = UniGene, 
           PubMed = PubMed, LocusLink = LocusLink )
               
        result.dir <- paste(dir, "Result", sep="/")
        dir.create(path=result.dir, showWarnings = FALSE)
        setwd(result.dir)
        normal.dir <- paste(result.dir, "Normalized_Data", sep="/")
        dir.create(path=normal.dir, showWarnings = FALSE)
        setwd(normal.dir)
        write.csv(Exprs, file="normaldata.csv", row.names=F)
        setwd(dir)
    }
    experimentData(normaldata)@preprocessing$normalMethod <- method
    return(normaldata)
}

