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
    control_code(50)
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
                                       
throwAwayTpts = 0
                
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

gamma = rep(0,nrow(S))

priorAlpha_gamma = 0.1
priorBeta_gamma = 1

beta = c(sim_results$true_fixed_beta, sim_results$true_time_varying_beta)
betaPriorPrecision = 1

N = matrix(data_list[["pop"]][,2], nrow = nrow(S), ncol = ncol(S), byrow = TRUE)



compMatDim = c(nrow(S), ncol(S))
outFileName = "./chainOutput_sim.txt"
# beta, rho,gamma, p_se, p_ei, p_ir,p_rs,S*,E*,I*,R*
logFileList = c(1, # beta
                1, # rho
                1, # gamma
                0, # p_se
                1, # p_ei
                1, # p_ir
                1, # p_rs 
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





iterationStride = 10
# S,E,R,beta,betaPrs,rho,gamma
sliceWidths = c(3,3,0.15,1e-4,1e-4,1e-4,1e-4)


# Check validity of dimensions

if (nrow(S) != nrow(E) || nrow(E) != nrow(I) || nrow(I) != nrow(R) || nrow(R) != nrow(S) || 
    ncol(S) != ncol(E) || ncol(E) != ncol(I) || ncol(I) != ncol(R) || ncol(R) != ncol(S) || 
    nrow(S) != length(gamma) || length(p_ei) != 1 || length(p_ir) != 1
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
                      sliceWidths,
                      reinfectionMode)

res$setRandomSeed(123123)
runSimulation = function(N, batchSize = 10)
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

tm = system.time(runSimulation(100))
print(tm)

