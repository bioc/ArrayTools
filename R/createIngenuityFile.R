`createIngenuityFile` <-
function(..., mydir=getwd(), eSet, filename="IngenuityFile"){

    dir <- mydir
    ##browser()
    all <- list(...)
    ##if (is.element(getAnnotation(all[[1]]), c("human.genechip", "mouse.genechip")) &&
    ##    is.null(annFile)) stop ("You need to provide the annotation file")

    result.dir <- paste(dir, "Result", sep="\\")
    dir.create(path=result.dir, showWarnings = FALSE)
    setwd(result.dir)
    path.dir <- paste(result.dir, "Pathway_Analysis", sep="\\")
    dir.create(path=path.dir, showWarnings = FALSE)
    setwd(path.dir)
    
    ing.dir <- paste(path.dir, "Ingenuity", sep="\\")
    dir.create(path=ing.dir, showWarnings = FALSE)   
    setwd(ing.dir)

    output.ing(all, eSet=eSet, filename=paste(ing.dir, filename, sep="\\"))
    setwd(dir)
}

