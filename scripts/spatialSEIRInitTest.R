library(spatialSEIR)
{if (file.exists("./simulation/SimulationObjects.robj"))
{
    cat("Using existing simulation data.\n")
    load("./simulation/SimulationObjects.robj")
}
else
{
    cat("Simulating data, this may take a moment.\n ")
    wd = getwd()
    setwd("./simulation")
    source("./simulateIowaData.R")
    control_code(20,10)
    load("./SimulationObjects.robj")
    setwd(wd)
}}


S = sim_results$S
E = sim_results$E
I = sim_results$I
R = sim_results$R

S_star = sim_results$S_star
E_star = sim_results$E_star
I_star = sim_results$I_star
R_star = sim_results$R_star


set.seed(123132);

X  = sim_results$X # X is returned by the simulation as an NxP matrix, where N 
                   # is the number of locations and P the number of predictors.
Z = sim_results$Z 
                                       
throwAwayTpts = 7*12
                
# Create basis for p_rs

X_betaPrs = sim_results$X_prs
X_betaPrs = X_betaPrs[(throwAwayTpts+1):nrow(S_star),]
X_betaPrsDim = dim(X_betaPrs)
betaPrs = sim_results$true_beta_pRS
betaPrsPriorPrecision = 1


S0 = S[throwAwayTpts + 1,]
E0 = E[throwAwayTpts + 1,]
I0 = I[throwAwayTpts + 1,]
R0 = R[throwAwayTpts + 1,]



R_star0 = R_star[throwAwayTpts,]

S_star = S_star[(throwAwayTpts+1):nrow(S_star),]
E_star = E_star[(throwAwayTpts+1):nrow(E_star),]
I_star = I_star[(throwAwayTpts+1):nrow(I_star),]
R_star = R_star[(throwAwayTpts+1):nrow(R_star),]
S = S[(throwAwayTpts+1):nrow(S),] 
E = E[(throwAwayTpts+1):nrow(E),] 
I = I[(throwAwayTpts+1):nrow(I),] 
R = R[(throwAwayTpts+1):nrow(R),] 

Z = Z[(1+throwAwayTpts*ncol(S)):nrow(Z),]

Z2 = c()

# Change order of Z
for (i in 1:ncol(S))
{
    Z2 = rbind(Z2, Z[seq(i, nrow(Z), ncol(S)),])
}
Z = Z2

xDim = dim(X)
zDim = dim(Z)

DM = as.numeric(data_list$dcm)

rho = 0.05

p_ei = 0.9
p_ir = 0.9

beta = c(sim_results$true_fixed_beta, sim_results$true_time_varying_beta)
betaPriorPrecision = 1

N = matrix(data_list[["pop"]][,2], nrow = nrow(S), ncol = ncol(S), byrow = TRUE)



compMatDim = c(nrow(S), ncol(S))
outFileName = "./chainOutput_sim.txt"


iterationStride = 100

# S,E,R,S0,I0,beta,betaPrs,rho
sliceWidths = c(0.2,  # S_star
                0.2,  # E_star
                0.2, # I_star
                0.2, # S0
                0.2, # I0
                0.001, # beta
                0.005, # betaPrs
                0.01# rho
                )




# Check validity of dimensions

if (nrow(S) != nrow(E) || nrow(E) != nrow(I) || nrow(I) != nrow(R) || nrow(R) != nrow(S) || 
    ncol(S) != ncol(E) || ncol(E) != ncol(I) || ncol(I) != ncol(R) || ncol(R) != ncol(S) || 
    length(p_ei) != 1 || length(p_ir) != 1
    || nrow(Z) != nrow(S)*ncol(S))
{
    stop("Invalid Starting Dimensions")
}

if (!all((S+E+I+R) == N) || any(S<0) || any(E<0) || any(I<0) ||
    any(R<0) || any(S_star<0) || any(E_star<0) || any(R_star<0))
{
    stop("Invalid Compartment Values")
}

verbose = FALSE
debug = FALSE

priorAlpha_pEI = 1000;
priorBeta_pEI = 100;
priorAlpha_pIR = 1000;
priorBeta_pIR = 100;

reinfectionMode = 1;
# Mode 1: estimate betaP_RS, S_star
# Mode 2: fix betaP_RS, estimate S_star
# Mode 3: no reinfection

# forget true values

true_beta = beta
true_beta_prs = betaPrs
true_rho = rho

beta = true_beta*0; beta[1] = 0.1
betaPrs = true_beta_prs*0; betaPrs[1] = 0.1
rho = 0.00001


steadyStateConstraintPrecision = 0.01

proposal = generateCompartmentProposal(I_star, N, S0, E0, I0)

res = spatialSEIRModel(compMatDim,
                      xDim,
                      zDim,
                      X_betaPrsDim,
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
                      X_betaPrs,
                      DM,
                      rho,
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
                      iterationStride,
                      steadyStateConstraintPrecision,
                      verbose,
                      debug, 
                      sliceWidths,
                      reinfectionMode)

# Use OpenCL:
res$samplingMode = 3
#res$oclPreferences = res$oclPreferences + 1 

res$setRandomSeed(123123)


Norder = order(N[1,])
itrPrint = function(x, wd=8)
{
    formatC(x, width = wd, format = "d", flag = "0")
}

makePlot = function(imgNo)
{
    png(filename = paste("./imgOut/", itrPrint(imgNo), ".png", sep =""), 
        width = 800, 
        height = 600) 
        
        plotTwoCompartments((I/N)[,Norder], 
                            (res$I/N)[,Norder], 
                            main1 = "True Infectious Proportion", 
                            main2 = "Fitted Infectious Proportion", zlim = c(0, min(max((I/N))*1.01, 1)))
    dev.off()
}

imgNo = 0
runSimulation = function(N, batchSize = 100, targetRatio = 0.25, targetWidth = 0.05, proportionChange = 0.01, printAR = FALSE)
{
    tryCatch({
        for (i in 1:(N/batchSize))
        {
            imgNo <<- imgNo + 1 
            makePlot(imgNo)
            res$simulate(batchSize)
            if (printAR)
            {
                res$printAcceptanceRates()
            }
            res$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
            # sleep to allow R to catch up and handle interrupts 
            Sys.sleep(0.001)
            cat(i*batchSize,"\n")

        }}, 
        interrupt = function(interrupt)
        {
            cat("Exiting...\n")
    })
}

#{
#targetRatio = 0.2
#targetWidth = 0.05
#proportionChange = 0.01
#for (j in 1:100)
#{
#    runSimulation(20, 1)
#    res$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
#}
#for (j in 1:10)
#{
#    runSimulation(1000, 1)
#    res$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
#}
#}

runSimulation(1000,10, printAR=TRUE)
runSimulation(50000,100, printAR=TRUE)
#runSimulation(10000,1, printAR=TRUE)
#runSimulation(10000,100, printAR=TRUE)
#runSimulation(10000000,1000, printAR=TRUE)

#


