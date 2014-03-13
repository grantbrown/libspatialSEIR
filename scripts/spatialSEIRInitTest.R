library(spatialSEIR)
{if (file.exists("./simulation/SimulationObjects.robj") &
    file.exists("./simulation/distmat.robj"))
{
    cat("Using existing simulation data.\n")
    load("./simulation/SimulationObjects.robj")
    load("./simulation/distmat.robj")
}
else
{
    cat("Simulating data, this may take a moment.\n ")
    wd = getwd()
    setwd("./simulation")
    source("./simulateIowaData.R")
    control_code(947)
    load("./SimulationObjects.robj")
    load("./distmat.robj")
    setwd(wd)
}}

X  = covariates$X # X is returned by the simulation as an NxP matrix, where N 
                  # is the number of locations and P the number of predictors.
Z = covariates$Z  # Z is returned by the simulation as an NxQxT1xT2 array, 
                  # where N is the number of locations, Q the number of
                  # time varying predictors, T1 the week number, and T2
                  # the year number. libSpatialSEIR doesn't use divided
                  # time indices like this, so we'll just read it into a 
                  # (N*T1*T2)xQ CovariateMatrix slot. 


# The compartmental "matrices" are returned by the simulation as 
# NxT1xT2 arrays. 
compMatDim = c(dim(sim_results$S[1]), prod(sim_results$S[2:3]))
xDim = dim(X)
zDim = c(prod(dim(Z)[c(1,3,4)]), dim(Z)[2])

S0 = -1
E0 = -1
I0 = -1
R0 = -1
Sstar0 = -1
Estar0 = -1
Istar0 = -1
Rstar0 = -1
Istar = sim_results$I_star



# Doesn't do anything interesting, just testing data transfer
res = spatialSEIRInit(compMatDim,
                      xDim,
                      zDim,
                      S0,
                      E0,
                      I0,
                      R0,
                      Sstar0,
                      Estar0,
                      Istar0,
                      Rstar0,
                      Istar,
                      X,
                      Z)

