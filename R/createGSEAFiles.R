`createGSEAFiles` <-
function(mydir=getwd(), eSet, catVar){
   
    dir <- mydir
    if (!is.element(catVar, colnames(pData(eSet))))     
        stop (paste(catVar, " is not found in the phenodata ", sep=""))

    result.dir <- paste(dir, "Result", sep="/")
    dir.create(path=result.dir, showWarnings = FALSE)
    setwd(result.dir)
    path.dir <- paste(result.dir, "Pathway_Analysis", sep="/")
    dir.create(path=path.dir, showWarnings = FALSE)
    setwd(path.dir)
    GSEA.dir <- paste(path.dir, "GSEA", sep="/")
    dir.create(path=GSEA.dir, showWarnings = FALSE)
    setwd(GSEA.dir)
    clsFileName <- paste(catVar, "phenotype", sep=".")
    output.cls(pData(eSet), catVar, paste(GSEA.dir, clsFileName , sep="/"))
    gctFileName <- paste(catVar, "probe", sep=".")
    output.gct(eSet, paste(GSEA.dir, gctFileName , sep="/"))
    setwd(dir)
}

