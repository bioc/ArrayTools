
geneFilter <- function (object, pct = 0.1, numChip = ceiling(ncol(exprs(object))*pct), 
    bg = 4, range = 0, iqrPct = 0, output =FALSE, mydir = getwd()) 
{

    if (bg > 0) {
        maxChip <- ncol(exprs(object))
        if (numChip > maxChip) {
            stop(paste("numChip must <= ", maxChip, sep = ""))
        }
        f <- genefilter::kOverA(numChip, bg)
        fun <- genefilter::filterfun(f)
        good1 <- genefilter::genefilter(object, fun)
    } else {
        good1 <- rep(TRUE, nrow(exprs(object)))
    }
    if (range > 0) {
        max <- apply(exprs(object), 1, max)
        min <- apply(exprs(object), 1, min)
        data_ratio <- max - min
        if (range > max(data_ratio)) {
            stop(paste("range must < ", max(data_ratio), sep = ""))
        }
        good2 <- data_ratio >= range
    } else {
        good2 <- rep(TRUE, nrow(exprs(object)))
    }
    if (iqrPct > 0) {
        iqr <- apply(exprs(object), 1, IQR)
        good3 <- iqr >= quantile(iqr, iqrPct)
    } else {
        good3 <- rep(TRUE, nrow(exprs(object)))
    }
    good <- good1 & good2 & good3
    print(paste("After Filtering, N = ", sum(good), sep = ""))
    filtered <- object[good,]
    if (output) {
        dir <- mydir
        packageName <- paste(annotation(object), ".db", sep="")
        library(packageName , character.only=TRUE)
        require("annaffy")      
        probeset_id = rownames(exprs(filtered))
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
        Exprs <- data.frame(probeset_id = probeset_id, exprs(filtered),
           Symbol = Symbol, Description = Description, Chromosome = Chromosome,
           GenBank = GenBank, Cytoband = Cytoband, UniGene = UniGene, 
           PubMed = PubMed, LocusLink = LocusLink )               
        result.dir <- paste(dir, "Result", sep="/")
        dir.create(path=result.dir, showWarnings = FALSE)
        setwd(result.dir)
        normal.dir <- paste(result.dir, "Filtered_Data", sep="/")
        dir.create(path=normal.dir, showWarnings = FALSE)
        setwd(normal.dir)
        write.csv(Exprs, file="filterdata.csv", row.names=F)
        setwd(dir)
    }
    experimentData(filtered)@other <- list(numChip = numChip, 
        bg =bg , range = range, iqrPct = iqrPct, nGene =dim(object)[1],
        nSample = dim(object)[2])
    filtered
}
