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

set.seed(123132);

X  = covariates$X # X is returned by the simulation as an NxP matrix, where N 
                  # is the number of locations and P the number of predictors.
Z_ar = covariates$Z  # Z is returned by the simulation as an NxQxT0xT2 array, 
                     # where N is the number of locations, Q the number of
                     # time varying predictors, T1 the week number, and T2
                     # the year number. libSpatialSEIR doesn't use divided
                     # time indices like this, so we'll just read it into a 
                     # (N*T1*T2)xQ CovariateMatrix slot. 
Z = Z_ar[,,1,1]

# Flatten the Z array 
for (id4 in 1:(dim(Z_ar)[4]))
{
    for (id3 in 1:(dim(Z_ar)[3]))
    {
        if (id4 != 1 || id3 != 1)
        {
            Z = rbind(Z, Z_ar[,,id3,id4])
        }
    }
}



# The compartmental "matrices" are returned by the simulation as 
# NxT1xT2 arrays. Define a function to convert these to 
# two dimensional Nx(T1xT2) matrices. 

flattenCompartment = function(compartment)
{
    output = compartment[,1,1]
    for (idx3 in 1:dim(compartment)[3])
    {
        for (idx2 in 1:dim(compartment)[2])
        {
            if (idx2 != 1 || idx3 != 1)
            {
                output = cbind(output, compartment[,idx2, idx3])
            }
        }
    }
    output
}

throwAwayTpts = 150

Sstar = flattenCompartment(sim_results$S_star)
Estar = flattenCompartment(sim_results$E_star)
Istar = flattenCompartment(sim_results$I_star)
Rstar = flattenCompartment(sim_results$R_star)

S = flattenCompartment(sim_results$S)
E = flattenCompartment(sim_results$E)
I = flattenCompartment(sim_results$I)
R = flattenCompartment(sim_results$R)


# Create basis for p_rs
Tvals = 1:ncol(S)
X_betaPrs = cbind(1,Tvals,sin(2*pi*(Tvals/52)), cos(2*pi*(Tvals/52)))
X_betaPrs = X_betaPrs[(throwAwayTpts+1):ncol(Sstar),]
X_betaPrsDim = dim(X_betaPrs)
betaPrs = c(3, rep(0, (ncol(X_betaPrs)-1)))
betaPrsPriorPrecision = 1




S0 = S[,throwAwayTpts]
E0 = E[,throwAwayTpts]
I0 = I[,throwAwayTpts]
R0 = R[,throwAwayTpts]
Sstar0 = Sstar[,throwAwayTpts]
Estar0 = Estar[,throwAwayTpts]
Istar0 = Istar[,throwAwayTpts]
Rstar0 = Rstar[,throwAwayTpts]

Sstar = Sstar[,(throwAwayTpts+1):ncol(Sstar)]
Estar = Estar[,(throwAwayTpts+1):ncol(Estar)]
Istar = Istar[,(throwAwayTpts+1):ncol(Istar)]
Rstar = Rstar[,(throwAwayTpts+1):ncol(Rstar)]
S = S[,(throwAwayTpts+1):ncol(S)] 
E = E[,(throwAwayTpts+1):ncol(E)] 
I = I[,(throwAwayTpts+1):ncol(I)] 
R = R[,(throwAwayTpts+1):ncol(R)] 

compMatDim = c(nrow(S), ncol(S))

Z = Z[(1+throwAwayTpts*nrow(S)):nrow(Z),]

xDim = dim(X)
zDim = dim(Z)

DM = as.numeric(data_list$dcm)

rho = 0.6

p_ei = 0.9
p_ir = 0.6
p_rs = rep(c(rep(0.1, 8), rep(0.1,4), rep(0.1,8), rep(0.4,4), rep(0.5, 3), 
         rep(0.9,3), rep(0.1, 4), rep(0.05, 18)), dim(sim_results$S)[3])
p_rs = p_rs[(throwAwayTpts+1):length(p_rs)]

gamma = rep(sim_results$gamma, dim(sim_results$S)[3])
gamma = gamma[(throwAwayTpts+1):length(gamma)]

priorAlpha_gamma = 0.1
priorBeta_gamma = 1

beta = c(covariates$true_fixed_beta, covariates$true_time_varying_beta)
N = matrix(data_list[["pop"]][,2], nrow = compMatDim[1], ncol = compMatDim[2])
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
                0) # Total pSE_j, j = 1...T




iterationStride = 10
# S,E,R,beta,betaPrs,rho,gamma
sliceWidths = c(3,3,3,1e-4,1e-4,1e-4,1e-4)


# Check validity of dimensions

if (nrow(S) != nrow(E) || nrow(E) != nrow(I) || nrow(I) != nrow(R) || nrow(R) != nrow(S) || 
    ncol(S) != ncol(E) || ncol(E) != ncol(I) || ncol(I) != ncol(R) || ncol(R) != ncol(S) ||
    ncol(S) != length(p_rs) || ncol(S) != length(gamma) || length(p_ei) != 1 || length(p_ir) != 1
    || nrow(Z) != nrow(S)*ncol(S))
{
    stop("Invalid Starting Dimensions")
}

if (!all((S+E+I+R) == N) || any(S<0) || any(E<0) || any(I<0) ||
    any(R<0) || any(Sstar<0) || any(Estar<0) || any(Rstar<0))
{
    stop("Invalid Compartment Values")
}

verbose = TRUE
debug = FALSE

priorAlpha_pEI = 1;
priorBeta_pEI = 1;
priorAlpha_pIR = 1;
priorBeta_pIR = 1;



res = spatialSEIRModel(compMatDim,
                      xDim,
                      zDim,
                      X_betaPrsDim,
                      S0,
                      E0,
                      I0,
                      R0,
                      Sstar0,
                      Estar0,
                      Istar0,
                      Rstar0,
                      Sstar,
                      Estar,
                      Istar,
                      Rstar,
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
                      betaPrs,
                      betaPrsPriorPrecision,
                      p_ei,
                      p_ir,
                      N,
                      outFileName, 
                      logFileList, 
                      iterationStride,
                      verbose,
                      debug, 
                      sliceWidths)


res$setRandomSeed(123123)
tryCatch({
    for (i in 1:100000)
    {
        res$simulate(1)
        Sys.sleep(0.001)
    }}, 
    interrupt = function(interrupt)
    {
        cat("Exiting...\n")
})


