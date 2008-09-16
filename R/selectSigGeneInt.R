`selectSigGeneInt` <-
function(object, pGroup = 0.05, fcGroup = 0, 
        pMain = 0.05, fcMain =0) {     
    sigResult.main <- selectSigGene(object[[2]], p.value = pMain,
        fc.value = fcMain)
    sigResult.group <- list()
    numStrata <- getLength(object) - 2
    for (i in 1:numStrata){
        sigResult.group[[i]] <- selectSigGene(object[[i+2]], p.value = pGroup,
            fc.value = fcGroup)
    }
    new("interactionResult", c(object[[1]], sigResult.main, sigResult.group))
    
}

