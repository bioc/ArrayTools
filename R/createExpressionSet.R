createExpressionSet <- function(pData, exprs, ...){
    exprs1 <- exprs[,grep(".CEL", colnames(exprs))]
    if (all(is.element(rownames(pData), colnames(exprs1))) &
        all(is.element(colnames(exprs1), rownames(pData)))) {
        pData1 <- pData[match(colnames(exprs1), rownames(pData)), ,drop=FALSE]
    } else stop ("The CEL file names in the expression data do not match 
        the ones in pheno data file/")
    phenoData <- new("AnnotatedDataFrame", data = pData1)
    eSet <- new("ExpressionSet", exprs = as.matrix(exprs1), phenoData = phenoData, ...)
    return(eSet)
}
