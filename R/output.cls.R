`output.cls` <-
function(target, variable, filename="phenotype"){       
    if(is.factor(target[,variable]))
        level.temp <- levels(target[,variable])
    else 
        level.temp <- levels(as.factor(target[,variable]))
    file <- as.vector(target[,variable])
    level <-level.temp[match(level.temp, unique(file))]
    filename1 <- paste(filename, ".cls", sep="")
    cat (c(length(file), length(level), 1), "\n", sep=" ", file=filename1)
    cat (c("#", level), "\n", sep=" ", file=filename1, append=T)
    cat (file, sep=" ", file=filename1, append=T)
}

