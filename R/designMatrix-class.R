
############ Design Matrix Class ####################

setClass("designMatrix", representation (design = "matrix", 
    target = "data.frame", covariates= "character", intIndex= "numeric"))
setMethod("getDesign", signature("designMatrix"), function(object) object@design)
setMethod("getTarget", signature("designMatrix"), function(object) object@target)
setMethod("getCovariates", signature("designMatrix"), function(object) object@covariates)
setMethod("getIntIndex", signature("designMatrix"), function(object) object@intIndex)
setMethod("show", signature("designMatrix"), function(object) {
    print(object@design) 
    invisible(object)
})
setMethod("initialize", signature("designMatrix"), function(.Object, ...,
    target, covariates, intIndex = 0) {

    if(!all(is.element(covariates, colnames(target))))
        stop ("The covariates are not matched with the names of the variables!")
    ##browser()
	
	for (i in 1:ncol(target)){
	    if (!is.factor(target[,i])) target[,i] <- as.factor(target[,i])
	}

    FORMULA <- createFormula(target, covariates, int = intIndex)
    design <- with(target, model.matrix(as.formula(FORMULA)))
    newcolnames <- colnames(design)
    covRearrange <- covariates[order(as.vector(sapply(covariates, nchar)), decreasing=T)]
    newcolnames1 <- newcolnames 

    for (i in 1:length(newcolnames1)){   
        if (length(grep(":", newcolnames1[i]))==0) {       
            for (j in 1:length(covRearrange )){
               if (length(grep("/", newcolnames1[i]))==0)            
                    newcolnames1[i] <- sub(covRearrange[j], paste(covRearrange [j], "/", sep = "") , 
                        newcolnames1[i])
            }
	   } else {
            var1 <- unlist(strsplit(newcolnames1[i], "\\:"))[1]
            var2 <- unlist(strsplit(newcolnames1[i], "\\:"))[2]
		for (j in 1:length(covRearrange )){
                if (length(grep("/", var1))==0)            
                    var1 <- sub(covRearrange[j], paste(covRearrange [j], "/", sep = "") , var1)
				if (length(grep("/", var2))==0) 
					var2 <- sub(covRearrange[j], paste(covRearrange [j], "/", sep = "") , var2)
            }
	      newcolnames1[i] <- paste(var1, var2, sep=":")
	  }		
    }

    colnames(design) <- newcolnames1
    callNextMethod(.Object, ...,  design = design, target = target, 
        covariates = covariates, intIndex= intIndex)
})

createFormula<-function (target, covariate, int)
{    
    factor.list <- list()
    index.cov <- formula.cov <- c()
    c.name<-colnames(target)
    if (!identical(int, 0)) {
        if (length(int) !=2) 
            stop ("You can only have two numbers")
        index.int1 <- which (colnames(target) == covariate[int[1]])
        index.int2 <- which (colnames(target) == covariate[int[2]])
        formula.int <- paste(c.name[index.int1],c.name[index.int2],sep="*")
        cov.not.int <- covariate[!is.element(covariate,covariate[int])]
        if (length(cov.not.int) > 0) {
            for (j in 1:length(cov.not.int)) {
		    index.cov[j] <- which (colnames(target) == cov.not.int[j])
		    formula.int <- paste(formula.int, c.name[index.cov[j]],sep="+")
		}
	  }
	  formula <- paste("~", formula.int, sep="")
    } else {
        index.cov[1] <- which (colnames(target) == covariate[1])
        formula.cov<-paste("~", c.name[index.cov[1]], sep="")
	  if (length(covariate) > 1) {
            for (j in 2:length(covariate)) {
                index.cov[j] <- which (colnames(target) == covariate[j])
                formula.cov<-paste(formula.cov, c.name[index.cov[j]], sep='+')
	        }
        }	
	  formula <- formula.cov
    }
    formula
}