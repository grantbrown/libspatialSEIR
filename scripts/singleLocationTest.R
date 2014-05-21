library(spatialSEIR)

set.seed(123123)
NYears = 20
TptPerYear = 12
MaxTpt = NYears*TptPerYear

ThrowAwayTpt = 180

X = matrix(1, ncol = 1)
Z = cbind(seq(1,NYears*TptPerYear), model.matrix(~as.factor(rep(1:12,NYears)))[,2:TptPerYear])

X_prs = cbind(1, 
              sin((1:MaxTpt)/TptPerYear*2*pi), 
              cos((1:MaxTpt)/TptPerYear*2*pi))


trueBetaSEFixed = c(-0.3)
trueBetaSEVarying = c(0.0001, 0.2, 0.3, 0.5, 0.6, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1)*3
trueGamma = rep(0.0, MaxTpt)  

eta_se = as.numeric((X %*% trueBetaSEFixed)) + (Z %*% trueBetaSEVarying)
p_se = numeric(MaxTpt)
p_ei = 0.9
p_ir = 0.9

trueBetaRS = c(2.5, -1, 0.5) 
eta_rs = X_prs %*% trueBetaRS
p_rs = exp(-eta_rs)

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

plotEpidemic = function()
{
    par(mfrow = c(3,1))
    plot(S[1,], type = "l", main = "S,E,I,R", ylim = c(0, max(S)))
    lines(E[1,], col = "red")
    lines(I[1,], col = "orange")
    lines(R[1,], col = "blue")

    plot(I[1,]/N, type = "l", col = "red", main = "Proportion Infected", ylim = c(0,1))
    plot(p_se, type = "l", col = "red", lty = 2, main = "Probability of Infection", ylim = c(0, max(p_se)))
}
plotEpidemic2 = function()
{
    par(mfrow = c(5,2))

    plot(S_star, type = "l", main = "S_star", ylim = c(0, max(max(S_star), max(res$S_star))))
    lines(res$S_star, col = "red")

    plot(S, type = "l", main = "S", ylim = c(0, max(max(S), max(res$S))))
    lines(res$S, col = "red")

    plot(E_star, type = "l", main = "E_star", ylim = c(0, max(max(E_star), max(res$E_star))))
    lines(res$E_star, col = "red")

    plot(E, type = "l", main = "E", ylim = c(0, max(max(E), max(res$E))))
    lines(res$E, col = "red")

    plot(I_star, type = "l", main = "I_star", ylim = c(0, max(max(I_star), max(res$I_star))))
    lines(res$I_star, col = "red")

    plot(I, type = "l", main = "I", ylim = c(0, max(max(I), max(res$I))))
    lines(res$I, col = "red")

    plot(R_star, type = "l", main = "R_star", ylim = c(0, max(max(R_star), max(res$R_star))))
    lines(res$R_star, col = "red")

    plot(R, type = "l", main = "R", ylim = c(0, max(max(R), max(res$R))))
    lines(res$R, col = "red")

    plot(p_se, type = "l", main = "p_se")
    lines(res$p_se, col = "red")

    plot(p_rs, type = "l", main = "p_rs") 
    lines(res$p_rs, col = "red")

}




#plotEpidemic()

# Format for libspatialSEIR

if (ThrowAwayTpt != 0)
{
    S0 = S[,ThrowAwayTpt]
    E0 = E[,ThrowAwayTpt]
    I0 = I[,ThrowAwayTpt]
    R0 = R[,ThrowAwayTpt]
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

# Not applicable parameters
DM = c(0)
rho = 0

p_ei = p_ei
p_ir = p_ir
p_rs = p_rs[(ThrowAwayTpt +1):length(p_rs)]
p_se = p_se[(ThrowAwayTpt +1):length(p_se)]
eta_se = eta_se[(ThrowAwayTpt +1):length(eta_se)]
gamma = trueGamma[(1+ThrowAwayTpt):length(trueGamma)]
beta = c(trueBetaSEFixed, trueBetaSEVarying) 
betaPrs = trueBetaRS
N = matrix(N, nrow = nrow(S), ncol = ncol(S))
outFileName = "./chainOutput_single.txt"
logFileList = c(1, # beta
                0, # rho
                0, # gamma
                0, # p_se
                1, # p_ei
                1, # p_ir
                0, # p_rs 
                0, # S*
                0, # E*
                0, # I*
                0, # R*
                1, # S total
                1, # E total
                1, # I total
                1, # R total
                1, # S_star total
                1, # E_star total
                1, # I_star total
                1, # R_star total
                1, # Average pSE
                1, # Average pRS
                0, # Total S_j, j = 1...T
                0, # Total E_j, j = 1...T
                0, # Total I_j, j = 1...T
                0, # Total R_j, j = 1...T
                0, # Total S_star_j, j = 1...T
                0, # Total E_star_j, j = 1...T
                0, # Total I_star_j, j = 1...T
                0, # Total R_star_j, j = 1...T
                0, # Total pSE_j, j = 1...T
                0, # Estimated R0, basic reproductive number
                0) # Estimated R0_j, j = 1...T 

iterationStride = 10000
# S,E,R,beta,betaPrs,rho,gamma
sliceWidths = c(5,5,5,1e-1,1e-1,1e-1,1e-1)

priorAlpha_gamma = 0.1
priorBeta_gamma = 1
priorAlpha_pEI = 1;
priorBeta_pEI = 1;
priorAlpha_pIR = 1;
priorBeta_pIR = 1;
betaPrsPriorPrecision = 0.1
betaPriorPrecision = 0.1

steadyStateConstraintPrecision = 0.001

verbose = FALSE 
debug = FALSE


# pretend not to know the true values of things
proposal = generateCompartmentProposal(I_star, N)
#beta = c(-1, rep(0, (length(beta)-1)))
#betaPrs = c(3, rep(0,(length(betaPrs)-1)))
#p_ei = 0.8
#p_ir = 0.8

res = spatialSEIRModel(compMatDim,
                      xDim,
                      zDim,
                      xPrsDim,
                      proposal$S0,
                      proposal$E0,
                      proposal$I0,
                      proposal$R0,
                      proposal$S_star,
                      proposal$E_star,
                      proposal$I_star,
                      proposal$R_star,
                      X,
                      Z,
                      X_prs,
                      DM,
                      rho,
                      gamma,
                      priorAlpha_gamma,
                      priorBeta_gamma,
                      priorAlpha_pEI,
                      priorBeta_pEI,
                      priorAlpha_pIR,
                      priorBeta_pIR,
                      beta,
                      betaPriorPrecision,
                      betaPrs,
                      betaPrsPriorPrecision,
                      p_ei,
                      p_ir,
                      N,
                      outFileName, 
                      logFileList, 
                      iterationStride,
                      steadyStateConstraintPrecision,
                      verbose,
                      debug, 
                      sliceWidths)


res$setRandomSeed(123123)

runSimulation = function(N, batchSize = 1000)
{
    tryCatch({
        for (i in 1:(N/batchSize))
        {
            res$simulate(batchSize)
            # sleep to allow R to catch up and handle interrupts 
            Sys.sleep(0.001)
            cat(i*batchSize,"\n")
        }}, 
        interrupt = function(interrupt)
        {
            cat("Exiting...\n")
    })
}

#runSimulation(10000)







