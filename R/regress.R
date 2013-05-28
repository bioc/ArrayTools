`regress` <-
function (object, contrast, method = c("limma", "regression", "permutation"), 
    adj = "none", permute.time = 1000) 
{
    method <- match.arg(method)
    fit <- lmFit(object, getDesign(contrast))
    fit2 <- contrasts.fit(fit, getContrast(contrast))
    if (method == "regression" | method == "permutation") {
        fit2$t <- fit2$coef / fit2$stdev.unscaled / fit2$sigma
        F.stat <- classifyTestsF(fit2, fstat.only = TRUE, df=fit2$df.residual)
        fit2$F <- as.vector(F.stat)
        df1 <- attr(F.stat, "df1")
        df2 <- attr(F.stat, "df2")
        fit2$F.p.value <- pf(fit2$F, df1, df2, lower.tail = FALSE)
        if (method == "permutation"){
            f <- fit2$F
            p <- matrix(NA, nrow = length(f), ncol = (permute.time - 1))
            for (i in 1:(permute.time - 1)) {
                p[, i] <- permute.1(object, design, contrast, f)
            }
            p.1 <- rep(1, length(f))
            p <- cbind(p, p.1)
            count <- apply(p, 1, sum)
            fit2$F.p.value <- count/permute.time
        }

    } else if (method == "limma") {
        fit2 <- eBayes(fit2)
    } 
    if (adj == "none")
        adj.P.Value <- p.adjust(fit2$F.p.value, method = "fdr")
    else 
        adj.P.Value <- p.adjust(fit2$F.p.value, method = adj)

    ##result<- list(as.vector(unlist(fit2$genes)), as.list(as.data.frame(fit2$coefficients)),
    ##    fit2$F, fit2$F.p.value, adj.P.Value, design, contrast, method, adj)
    FC <- as.list(as.data.frame(fit2$coefficients))
    Log2Ratio.name <- character(ncol(getContrast(contrast)))
    for (i in 1:ncol(getContrast(contrast))) {
        Log2Ratio.name[i] <- paste("Log2Ratio", i, sep = ".")
    }
    names(FC) <- Log2Ratio.name
    
    result <- new("regressResult", ID = as.vector(unlist(rownames(fit2))), 
        foldChange = FC, FValue = fit2$F, pValue = fit2$F.p.value, 
        adjPVal = adj.P.Value, contrast =contrast,
        regressionMethod = method, adjustment = adj, annotation=object@annotation,
		normalizationMethod = object@experimentData@preprocessing, 
		filterMethod = object@experimentData@other)
    return(result)
}

