library(spatialSEIR)
set.seed(12313)

dataSubsetStartDate = "2008-01-01"
dataSubsetEndDate = "2009-06-01"


####################################
# Part 1: Read in and Process Data #
####################################

# Code below copied pretty directly from some past work. 
# Please excuse the mess. 

# Influenza Data - Via Google Flu Trends
# Originally used in my preceptorship project. 
fluData = read.csv("./scriptData/fluTrends/fluDataProcessed.csv")
fluDates = as.Date(as.character(fluData[,1]), format = "%m/%d/%Y")

fluData[,1] = fluDates
keepDates = (fluDates >= as.Date(dataSubsetStartDate, format = "%Y-%m-%d") & 
             fluDates < as.Date(dataSubsetEndDate, format = "%Y-%m-%d"))
fluData = fluData[keepDates,]
fluDates = fluDates[keepDates]

#fluMonthDate for matching to temperature Data
fluDateParts = strsplit(as.character(fluDates), "-")
fluMonthDate = lapply(fluDateParts, function(x){as.Date(paste(x[1], x[2], "01",
                                            sep = "-"), format = "%Y-%m-%d")})
fluMonthDate = as.Date(as.numeric(fluMonthDate), origin = "1970-01-01")
monthVal = c()
for (mth in strsplit(as.character(fluMonthDate), "-"))
{
    monthVal = c(monthVal, as.numeric(mth[2]))
}


Y = as.matrix(fluData[,-1])
Y = Y[,order(colnames(Y))]
Y = t(Y)

yMatchID = paste(rownames(Y), as.character(fluMonthDate), sep = "-")


# Population data via the US Census Bureau 
fluPopulation = read.csv("./scriptData/fluTrends/fluPopulation.csv")
rownames(fluPopulation) = fluPopulation$State
fluPopulation = fluPopulation[order(rownames(fluPopulation)),]

# Temperature Data
temperature = read.csv("./scriptData/fluTrends/temperatureData.csv")
temperature$Date = as.Date(paste(as.character(temperature$Year),
                                 as.character(temperature$Month),
                                 rep("01", nrow(temperature)),
                                 sep = "-"), format = "%Y-%m-%d")
tempKeepDates = (temperature$Date >= as.Date(dataSubsetStartDate, format = "%Y-%m-%d") & 
                 temperature$Date < as.Date(dataSubsetEndDate, format = "%Y-%m-%d"))
temperature = temperature[tempKeepDates,]


temperature = data.frame(State = temperature$State, Date = temperature$Date, Temp = temperature$Temp)
temperature$uqId = paste(as.character(temperature$State), as.character(temperature$Date), sep = "-")


# Read in neighborhood information
neighborhood = as.matrix(read.csv("./scriptData/fluTrends/neighborhood.csv"))
rownames(neighborhood) = colnames(neighborhood)
neighborhood = neighborhood[order(colnames(neighborhood)), order(colnames(neighborhood))]


#####################################
# Part 2: Format for libSpatialSEIR #
#####################################

# Prepare Covariates
X  = cbind(1, (fluPopulation$Prop-mean(fluPopulation$Prop)))  
Z_ar = array(0, dim=c(nrow(Y),ncol(Y),1)) # Just a matrix in this case
# Multiple observations per monthly temperature adjustment
idx = c()
idxMult = cbind(1:length(unique(fluMonthDate)), table(fluMonthDate))
for (i in 1:nrow(idxMult))
{
    idx = c(idx, rep(idxMult[i,1], idxMult[i,2]))
}
states = unique(temperature$State)
for (i in 1:length(states))
{
    stateSub = temperature[temperature$State == states[i],]
    Z_ar[i,,1] = (stateSub$Temp - mean(temperature$Temp))[idx]
}

Z = matrix(Z_ar[,1,], ncol = 1)

# Flatten the Z array 
for (idx in 2:(dim(Z_ar)[2]))
{
    Z = rbind(Z, matrix(Z_ar[,idx,], ncol = 1)) 
}
Z = Z - mean(Z)


# Create N Matrix (not time varying at this point)
N = matrix(fluPopulation$Pop, nrow = nrow(Y), ncol = ncol(Y))

# Guess Initial Compartments. 

I_star = floor(sqrt(Y))

S0 = floor(0.95*N[,1]) 
E0 = I_star[,1]
I0 = floor(sqrt(Y[,1])/2)
R0 = N[,1] - floor(S0 + E0 + I0)

S = E = I = R = S_star = E_star = R_star = I_star*0

for (tpt in 1:ncol(R))
{
    if (tpt == 1)
    {
        S[,1] = S0
        E[,1] = E0
        I[,1] = I0
        R[,1] = R0
        S_star[,1] = rbinom(rep(1, length(S0)), R[,1], 0.05)
        E_star[,1] = rbinom(rep(1, length(S0)), I_star[,2], 1)
        # I_star fixed
        R_star[,1] = rbinom(rep(1, length(S0)), I[,1], 0.8)
    }
    else
    {
        S[,tpt] = S[,tpt-1] + S_star[,tpt-1] - E_star[,tpt-1]
        E[,tpt] = E[,tpt-1] + E_star[,tpt-1] - I_star[,tpt-1]
        I[,tpt] = I[,tpt-1] + I_star[,tpt-1] - R_star[,tpt-1]
        R[,tpt] = R[,tpt-1] + R_star[,tpt-1] - S_star[,tpt-1]

        S_star[,tpt] = rbinom(rep(1,length(S0)), R[,tpt], 0.05)
        if (tpt != ncol(R))
        {
            E_star[,tpt] = rbinom(rep(1,length(S0)), I_star[,tpt+1], 1)
            #E_star[,tpt] = I_star[,tpt+1]
        }
        else
        {
            E_star[,tpt] = E_star[,tpt-1]
        }
        # I_star fixed
        R_star[,tpt] = rbinom(rep(1,length(S0)), I[,tpt], 0.8)
    }
}


# Transpose everything for libspatialSEIR

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

N = t(N)



X_betaPrs = model.matrix(~as.factor(monthVal))[,1:11]
betaPrs = rep(0,ncol(X_betaPrs)) 
betaPrs[1] = 2.5
betaPrsPriorPrecision = 1
xDim = dim(X)
zDim = dim(Z)
X_betaPrsDim = dim(X_betaPrs)
compMatDim = c(nrow(S), ncol(S))

DM = neighborhood

rho = 0.1

p_ei = 0.99
p_ir = 0.8
p_rs = rep(0.1, nrow(S))
gamma = rep(0.1, nrow(S))
priorAlpha_gamma = 0.1
priorBeta_gamma = 1


#beta = c(-2,0.1,0.1)
beta = c(0,0,0)
betaPriorPrecision = 1

priorAlpha_pEI = 1000
priorBeta_pEI = 100
priorAlpha_pIR = 1000
priorBeta_pIR = 100

outFileName = "./chainOutput.txt"
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
                0, # Estimated R0, basic reproductive number. 
                0) # Estimated R0_j, j=1...T




iterationStride = 10

# S,E,R,beta,betaPrs,rho,gamma
sliceWidths = c(0.4,0.3,0.05,0.25,0.25,0.015,0.01)

if (!all((S+E+I+R) == N) || any(S<0) || any(E<0) || any(I<0) ||
    any(R<0) || any(S_star<0) || any(E_star<0) || any(R_star<0))
{
    stop("Invalid Compartment Values")
}
# Output Options
verbose = FALSE
debug = FALSE

reinfectionMode = 1
# Mode 1: estimate betaP_RS, S_star
# Mode 2: fix betaP_RS, estimate S_star
# Mode 3+: No reinfection

steadyStateConstraintPrecision = 0.01

res = spatialSEIRModel(compMatDim,
                      xDim,
                      zDim,
                      X_betaPrsDim,
                      S0,
                      E0,
                      I0,
                      R0,
                      S_star,
                      E_star,
                      I_star,
                      R_star,
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

runSimulation = function(N, batchSize = 100)
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

runSimulation(1000,10)
res$printAcceptanceRates()




