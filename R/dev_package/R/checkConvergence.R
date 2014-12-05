checkConvergence = function(fileName1,fileName2,fileName3,maxVal=2,useUpper=FALSE,verbose=TRUE)
{
    dat1 = read.csv(fileName1)
    dat2 = read.csv(fileName2)
    dat3 = read.csv(fileName3)
    maxIdx = min(nrow(dat1), nrow(dat2), nrow(dat3))
    dat1_sub = dat1[1:maxIdx,]
    dat2_sub = dat2[1:maxIdx,]
    dat3_sub = dat3[1:maxIdx,]


    
    iteration = dat1$Iteration

    varNames = names(dat1)
    lastVar = which(varNames == "gamma_ir") 
    dat1_mc = as.mcmc(dat1_sub[,1:lastVar])
    dat2_mc = as.mcmc(dat2_sub[,1:lastVar])
    dat3_mc = as.mcmc(dat3_sub[,1:lastVar])
    mcl = mcmc.list(dat1_mc,dat2_mc,dat3_mc)   

    

    diag = tryCatch({gelman.diag(mcl, multivariate=FALSE)}
        , warning=function(w){return(FALSE)}
        , error = function(e){return(FALSE)}
    )
    if (class(diag) == "logical" && !diag)
    {
        cat("WARNING: Convergence Check Failed. Continuing anyway...\n")
        return(FALSE)
    }

    if (useUpper)
    {
        criterion = max(diag[[1]][,2])
        maxParam = rownames(diag[[1]])[which.max(diag[[1]][,2])]
    }
    else
    {
        criterion = max(diag[[1]][,1])
        maxParam = rownames(diag[[1]])[which.max(diag[[1]][,1])]
    }
    if (verbose)
    {
        cat(paste("Convergence Criterion: ", round(criterion, 2), 
                  ", param: ", maxParam, "\n", sep = ""))
    }
    criterion < maxVal
}


