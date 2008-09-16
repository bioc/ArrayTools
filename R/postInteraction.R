`postInteraction` <-
function (eSet, regressObject, mainVar, compare1, compare2, 
    method=regressionMethod(regressObject), adj = adjustment(regressObject)) 
{
    
    if (class(regressObject) != "regressResult")
        stop ("The regressObject needs to be a 'regressResult' class")
    if (class(eSet) != "ExpressionSet")
        stop ("The eSet needs to be a ExpressionSet")

    covariates <- getCovariates(getContrast(regressObject))
    mainVarIndex <- which(mainVar == covariates)

    if (length(mainVarIndex) == 0)
        stop ("The main.var is not correct - it's not in the covariate list")
    
    strataVar <- covariates[getIntIndex(getContrast(regressObject))][-mainVarIndex]
    target <- pData(eSet)
   
    contrast.sep <- list()
    result.list <- list()
    u.1 <- as.vector(t(unique(target[strataVar])))
    design <- new("designMatrix", target = target, covariates=covariates,
        intIndex = getIntIndex(getContrast(regressObject)))

    for (j in 1:length(u.1)) {               
        contrast.sep[[j]] <- new("contrastMatrix", design.matrix = design, 
            compare1=compare1, compare2=compare2, level = u.1[j])
        result.list[[j]] <- regress(eSet[getIndex(regressObject),], contrast.sep[[j]], 
            method = method, adj = adj)
    }
    
    adjustVar <- covariates[-getIntIndex(getContrast(regressObject))]
    mainVar1 <- c(mainVar, adjustVar)
    design.main <- new("designMatrix", target=target, covariates=mainVar1)
    contrast.main <- new("contrastMatrix", design.matrix = design.main,
        compare1 = compare1, compare2 = compare2)
    result.main <- regress(eSet[!getIndex(regressObject),], contrast.main, 
            method = method, adj = adj)
    ##return (list(result.main, result.list))
    ##browser()
    new("interactionResult", c(regressObject, result.main, result.list))

}

