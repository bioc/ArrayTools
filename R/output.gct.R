`output.gct` <-
function(normal, filename="probe"){		
    NAME<-Description<-rownames(exprs(normal))
    file<-cbind(NAME, Description, exprs(normal))
    filename1 <- paste(filename, ".gct", sep="")
    cat("#1.2", "\n", sep="\t", file=filename1)
    cat(nrow(normal), ncol(normal), "\n", sep="\t", file=filename1, append=T)
    suppressWarnings(write.table(file, file=filename1, row.names=FALSE, 
        quote = FALSE, sep="\t", append=T))	
}

