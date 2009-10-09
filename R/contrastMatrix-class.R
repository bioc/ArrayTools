
################ Contrast Matrix Class #################
setClass("contrastMatrix", contains = "designMatrix", 
    representation (contrast = "matrix", 
    compare1 = "character", compare2 = "character",
    level = "character", interaction = "logical"))


setMethod("getContrast", signature("contrastMatrix"), function(object) object@contrast)
setMethod("getCompare1", signature("contrastMatrix"), function(object) object@compare1)
setMethod("getCompare2", signature("contrastMatrix"), function(object) object@compare2)
setMethod("getLevel", signature("contrastMatrix"), function(object) object@level)
setMethod("getInteraction", signature("contrastMatrix"), function(object) object@interaction)
setMethod("show", signature("contrastMatrix"), function(object) {
    print(object@contrast) 
    invisible(object)
})


setMethod("initialize", signature("contrastMatrix"), function(.Object, ..., 
    design.matrix = new("designMatrix"), compare1 = NA, compare2 = NA, level = NA, 
    interaction = FALSE) {
    
    param <- colnames(getDesign(design.matrix))
    if (interaction) {
        inter.eq <- array(rep(0, length(param) * length(param)), 
            dim = c(length(param), length(param)))
        inter.index <- rep(FALSE, length(param))
        for (j in 2:length(param)) {
            if (length(grep(":", param[j])) != 0) {
                inter.eq[j, j] <- 1
                inter.index[j] <- TRUE
            }
        }
        eq <- cbind(inter.eq[, inter.index])
    } else {
        eq1 <- eq2 <- c(1, rep(0, length(param)-1))
        for (i in 2:length(param)) {

            if (length(grep(":", param[i])) == 0) {
                value <- unlist(strsplit(param[i], "\\/"))[2]
                if (compare1 == value) 
                  eq1[i] <- 1
                if (compare2 == value) 
                  eq2[i] <- 1
            } else {
                var.1 <- unlist(strsplit(param[i], "\\:"))[1]
                var.2 <- unlist(strsplit(param[i], "\\:"))[2]
                v1<-unlist(strsplit(var.1, "\\/"))[2]
                v2<-unlist(strsplit(var.2, "\\/"))[2]              
                eq1[i] <- ifelse (compare1 == v1 | level == v1, 1, 0)*
		        ifelse (compare1 == v2 | level == v2, 1, 0)
                eq2[i] <- ifelse (compare2 == v1 | level == v1, 1, 0)*
			  ifelse (compare2 == v2 | level == v2, 1, 0)
            }
        }
        eq <- cbind(eq1 - eq2)
    }

    .Object <- callNextMethod(.Object, ..., 
        target = getTarget(design.matrix), covariates = getCovariates(design.matrix),
            intIndex = getIntIndex(design.matrix))
    .Object@contrast <- eq
    .Object@compare1 <- as.character(compare1)
    .Object@compare2 <- as.character(compare2)
    .Object@level <- as.character(level)
    .Object@interaction <- interaction
    .Object   
    
})