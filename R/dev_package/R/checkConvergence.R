checkConvergence = function(fileName1,fileName2,fileName3,maxVal=2,useUpper=FALSE,verbose=TRUE)
{
    dat1 = read.csv(fileName1)
    dat2 = read.csv(fileName2)
    dat3 = read.csv(fileName3)
    
    iteration = dat1$Iteration

    varNames = names(dat1)
    lastVar = which(varNames == "gamma_ir") 
    dat1 = as.mcmc(dat1[,1:lastVar])
    dat2 = as.mcmc(dat2[,1:lastVar])
    dat3 = as.mcmc(dat3[,1:lastVar])
    mcl = mcmc.list(dat1,dat2,dat3)
    diag = gelman.diag(mcl, multivariate=FALSE)
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


