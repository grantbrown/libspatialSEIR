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
    control_code(947)
    load("./SimulationObjects.robj")
    setwd(wd)
}}

X  = covariates$X # X is returned by the simulation as an NxP matrix, where N 
                  # is the number of locations and P the number of predictors.
Z_ar = covariates$Z  # Z is returned by the simulation as an NxQxT1xT2 array, 
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

compMatDim = c(dim(sim_results$S)[1], prod(dim(sim_results$S)[2:3]))
xDim = dim(X)
zDim = dim(Z)

S0 = sim_results$S0
E0 = sim_results$E0
I0 = sim_results$I0
R0 = sim_results$R0
Sstar0 = sim_results$S_star0
Estar0 = sim_results$E_star0
Istar0 = sim_results$I_star0
Rstar0 = sim_results$R_star0
Sstar = flattenCompartment(sim_results$S_star)
Estar = flattenCompartment(sim_results$E_star)
Istar = flattenCompartment(sim_results$I_star)
Rstar = flattenCompartment(sim_results$R_star)
DM = as.numeric(data_list$dcm)

rho = 0.6

p_ei = 0.8
p_ir = 0.6
p_rs = rep(c(rep(0.1, 8), rep(0.1,4), rep(0.1,8), rep(0.4,4), rep(0.5, 3), 
         rep(0.9,3), rep(0.1, 4), rep(0.05, 18)), dim(sim_results$S)[3])

beta = c(covariates$true_fixed_beta, covariates$true_time_varying_beta)
N = data_list[["pop"]][,2]

res = spatialSEIRInit(compMatDim,xDim,
                      zDim,S0,
                      E0,I0,
                      R0,Sstar0,
                      Estar0,Istar0,
                      Rstar0,Sstar,
                      Estar,Istar,
                      Rstar,X,
                      Z,DM,
                      rho,beta,
                      p_ei,p_ir,
                      p_rs,N)


