library(spatialSEIR)

set.seed(123123)
NYears = 3
TptPerYear = 12 
MaxTpt = NYears*TptPerYear
ThrowAwayTpt = 0


plotEpidemic = function(simResults,res)
{
    par(mfrow = c(5,2))

    plot(simResults$S_star, type = "l", main = "S_star", ylim = c(0, max(max(simResults$S_star), max(res$S_star))))
    lines(res$S_star, col = "red")

    plot(simResults$S, type = "l", main = "S", ylim = c(0, max(max(simResults$S), max(res$S))))
    lines(res$S, col = "red")

    plot(simResults$E_star, type = "l", main = "E_star", ylim = c(0, max(max(simResults$E_star), max(res$E_star))))
    lines(res$E_star, col = "red")

    plot(simResults$E, type = "l", main = "E", ylim = c(0, max(max(simResults$E), max(res$E))))
    lines(res$E, col = "red")

    plot(simResults$I_star, type = "l", main = "I_star", ylim = c(0, max(max(simResults$I_star), max(res$I_star))))
    lines(res$I_star, col = "red")

    plot(simResults$I, type = "l", main = "I", ylim = c(0, max(max(simResults$I), max(res$I))))
    lines(res$I, col = "red")

    plot(simResults$R_star, type = "l", main = "R_star", ylim = c(0, max(max(simResults$R_star), max(res$R_star))))
    lines(res$R_star, col = "red")

    plot(simResults$R, type = "l", main = "R", ylim = c(0, max(max(simResults$R), max(res$R))))
    lines(res$R, col = "red")

    plot(simResults$p_se, type = "l", main = "p_se")
    lines(res$p_se, col = "red")

    plot(simResults$p_rs, type = "l", main = "p_rs", ylim = c(0,1)) 
    lines(res$p_rs, col = "red")

}


runSimulation = function()
{
    X = matrix(1, ncol = 1)
    Z = cbind(seq(1,NYears*TptPerYear), sin(seq(1,NYears*TptPerYear)/TptPerYear*2*pi))
    X_prs = cbind(1, 
                  sin((1:MaxTpt)/TptPerYear*2*pi), 
                  cos((6+1:MaxTpt)/TptPerYear*2*pi))

    trueBetaSEFixed = c(0.5)
    trueBetaSEVarying = c(0.002, .5)
    trueGamma = rep(0.0, MaxTpt)  

    eta_se = as.numeric((X %*% trueBetaSEFixed)) + (Z %*% trueBetaSEVarying)
    p_se = numeric(MaxTpt)
    gamma_ei = 2.3
    gamma_ir = 2.3

    p_ei = 1-exp(-gamma_ei)
    p_ir = 1-exp(-gamma_ir)

    trueBetaRS = c(-2.5, 1, -0.25) 
    eta_rs = X_prs %*% trueBetaRS
    p_rs = 1-exp(-exp(eta_rs))

    N = 100000
    E0 = 0
    I0 = floor(0.001*N)
    R0 = floor(0.001*N) 
    S0 = N-E0-I0-R0

    S = matrix(0, nrow = 1, ncol = MaxTpt)
    E = matrix(0, nrow = 1, ncol = MaxTpt)
    I = matrix(0, nrow = 1, ncol = MaxTpt)
    R = matrix(0, nrow = 1, ncol = MaxTpt)
    S_star = matrix(0, nrow = 1, ncol = MaxTpt)
    E_star = matrix(0, nrow = 1, ncol = MaxTpt)
    I_star = matrix(0, nrow = 1, ncol = MaxTpt)
    R_star = matrix(0, nrow = 1, ncol = MaxTpt)



    for (i in 1:MaxTpt)
    {
        if (i == 1)
        {
            S[i] = S0
            E[i] = E0
            I[i] = I0
            R[i] = R0

            etaVal = -(I[i]/N)*exp(eta_se[i]) - trueGamma[i]
            p_se[i] = 1-exp(etaVal)

            S_star[i] = rbinom(1, R[i], p_rs[i])
            E_star[i] = rbinom(1, S[i], p_se[i])
            I_star[i] = rbinom(1, E[i], p_ei)
            R_star[i] = rbinom(1, I[i], p_ir)
        }
        else
        {
            S[i] = S[i-1] + S_star[i-1] - E_star[i-1]
            E[i] = E[i-1] + E_star[i-1] - I_star[i-1]
            I[i] = I[i-1] + I_star[i-1] - R_star[i-1]
            R[i] = R[i-1] + R_star[i-1] - S_star[i-1]

            etaVal = -(I[i]/N)*exp(eta_se[i]) - trueGamma[i]
            p_se[i] = 1-exp(etaVal)

            S_star[i] = rbinom(1, R[i], p_rs[i])
            E_star[i] = rbinom(1, S[i], p_se[i])
            I_star[i] = rbinom(1, E[i], p_ei)
            R_star[i] = rbinom(1, I[i], p_ir)
        }
    }

    if (ThrowAwayTpt != 0)
        {
        S0 = S[,ThrowAwayTpt+1]
        E0 = E[,ThrowAwayTpt+1]
        I0 = I[,ThrowAwayTpt+1]
        R0 = R[,ThrowAwayTpt+1]
    }

    S_star = S_star[,(ThrowAwayTpt + 1):ncol(S_star), drop = FALSE]
    E_star = E_star[,(ThrowAwayTpt + 1):ncol(E_star), drop = FALSE]
    I_star = I_star[,(ThrowAwayTpt + 1):ncol(I_star), drop = FALSE]
    R_star = R_star[,(ThrowAwayTpt + 1):ncol(R_star), drop = FALSE]

    S = S[,(ThrowAwayTpt + 1):ncol(S), drop = FALSE]
    E = E[,(ThrowAwayTpt + 1):ncol(E), drop = FALSE]
    I = I[,(ThrowAwayTpt + 1):ncol(I), drop = FALSE]
    R = R[,(ThrowAwayTpt + 1):ncol(R), drop = FALSE]

    Z = Z[(1+ThrowAwayTpt*nrow(S)):nrow(Z),]
    X_prs = X_prs[(ThrowAwayTpt + 1):nrow(X_prs),]

# Transpose Everything to have TXP

    S0 = t(S0)
    E0 = t(E0)
    I0 = t(I0)
    R0 = t(R0)

    S_star = t(S_star)
    E_star = t(E_star)
    I_star = t(I_star)
    R_star = t(R_star)

    S = t(S)
    E = t(E)
    I = t(I)
    R = t(R)

    xDim = dim(X)
    zDim = dim(Z)
    xPrsDim = dim(X_prs)
    compMatDim = c(nrow(S), ncol(S))

    p_rs = p_rs[(ThrowAwayTpt +1):length(p_rs)]
    p_se = p_se[(ThrowAwayTpt +1):length(p_se)]
    eta_se = eta_se[(ThrowAwayTpt +1):length(eta_se)]
    beta = c(trueBetaSEFixed, trueBetaSEVarying) 
    betaPrs = trueBetaRS
    N = matrix(N, nrow = nrow(S), ncol = ncol(S))

    return(list("S_star"=S_star, "E_star" = E_star, "I_star" = I_star, "R_star" = R_star,
                "S"=S, "E"=E, "I"=I, "R"=R, "S0"=S0,"E0"=E0,"I0"=I0,"R0"=R0,"X"=X, "Z"=Z,"X_prs"=X_prs,
                "p_se"=p_se,"p_rs"=p_rs, "beta" = beta, "betaPrs"=betaPrs, "N"=N))
}

simResults = runSimulation()

DataModel = buildDataModel(simResults$I_star, type = "overdispersion", params = c(1,1))
ExposureModel = buildExposureModel(simResults$X, simResults$Z, 
                                   beta = c(2, rep(0, ((length(simResults$beta))-1))), betaPriorPrecision = 0.1)
ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = simResults$X_prs, 
                                         betaPrs = -c(4, rep(0,(length(simResults$betaPrs)-1))), 
                                         priorPrecision = 0.1)
SamplingControl = buildSamplingControl(iterationStride=1000,
                                       sliceWidths = c(0.26,  # S_star
                                                       0.1,  # E_star
                                                       0.15, # I_star
                                                       0.22, # S0
                                                       0.24, # I0
                                                       0.8, # beta
                                                       0.2, # betaPrs
                                                       0.015, # rho
                                                       0.01, # gamma_ei
                                                       0.01, # gamma_ir
                                                       0.01 # phi
                                                      ))
DistanceModel = buildDistanceModel(list(matrix(0)))
TransitionPriors = buildTransitionPriorsManually(2300,1000,2300,1000)
InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, S0 = simResults$N[1]-100, 
                                           I0 = 100, E0 = 0)

outFileName = "./chainOutput_single.txt"

res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,TransitionPriors,
                     InitContainer,SamplingControl)

res$setRandomSeed(123123)
itrPrint = function(x, wd=8)
{
    formatC(x, width = wd, format = "d", flag = "0")
}


imgNo = 0;
runSimulation = function(N, batchSize = 100, targetRatio = 0.15, targetWidth = 0.05, proportionChange = 0.1, printAR = FALSE)
{
    tryCatch({
        for (i in 1:(N/batchSize))
        {
            imgNo <<- imgNo + 1 
            res$simulate(batchSize)
            if (printAR)
            {
                res$printAcceptanceRates()
            }
            res$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
            # sleep to allow R to catch up and handle interrupts 
            Sys.sleep(0.001)
            png(filename = paste("./imgOut/single/", itrPrint(imgNo), ".png", sep =""), width = 600, height = 1200) 
                plotEpidemic(simResults,res)
            dev.off()
            cat(i*batchSize,"\n")
        }}, 
        interrupt = function(interrupt)
        {
            cat("Exiting...\n")
    })
}

res$compartmentSamplingMode = 14
runSimulation(20000,100, printAR = FALSE, targetRatio = 0.2)
runSimulation(10000000,10000, printAR = TRUE, targetRatio = 0.2, targetWidth = 0.05)





