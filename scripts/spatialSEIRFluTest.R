library(spatialSEIR)

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
keepDates = (fluDates >= as.Date("2008-01-01", format = "%Y-%m-%d") & 
             fluDates < as.Date("2011-06-01", format = "%Y-%m-%d"))
fluData = fluData[keepDates,]
fluDates = fluDates[keepDates]

#fluMonthDate for matching to temperature Data
fluDateParts = strsplit(as.character(fluDates), "-")
fluMonthDate = lapply(fluDateParts, function(x){as.Date(paste(x[1], x[2], "01",
                                            sep = "-"), format = "%Y-%m-%d")})
fluMonthDate = as.Date(as.numeric(fluMonthDate), origin = "1970-01-01")



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
tempKeepDates = (temperature$Date >= as.Date("2008-01-01", format = "%Y-%m-%d") & 
                 temperature$Date < as.Date("2011-06-01", format = "%Y-%m-%d"))
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
X  = cbind(1, fluPopulation$Prop)  
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

Z = as.numeric(Z_ar)


# Create N Matrix (not time varying at this point)
N = matrix(fluPopulation$Pop, nrow = nrow(Y), ncol = ncol(Y))

# Guess Initial Compartments. 
E_star = Y[,c(2:(ncol(Y)-1),1)]
I_star = Y

S = E = I = R = S_star = R_star = I_star*0

for (loc in 1:nrow(R))
{
    for (tpt in 1:ncol(R))
    {
        if (tpt == 1)
        {
            S[loc,tpt] = min(floor(0.95*N[loc,tpt]), N[loc,tpt]-I_star[loc,tpt])
            E[loc,tpt] = E_star[loc, tpt] + Y[loc,tpt] 
            I[loc, tpt] = I_star[loc,tpt] 
            R[loc,tpt] = N[loc,tpt] - S[loc,tpt] - E[loc, tpt] - I[loc,tpt]              
        }
        S_star[loc,tpt] = rbinom(1,R[loc,tpt], 0.1)
        R_star[loc,tpt] = rbinom(1,I[loc,tpt], 0.8)
        if (tpt != ncol(R))
        {
            S[loc, tpt+1] = S[loc, tpt] + S_star[loc,tpt] - E_star[loc,tpt]
            E[loc, tpt+1] = E[loc, tpt] + E_star[loc,tpt] - I_star[loc,tpt]
            I[loc, tpt+1] = I[loc, tpt] + I_star[loc,tpt] - R_star[loc,tpt]
            R[loc, tpt+1] = R[loc, tpt] + R_star[loc,tpt] - S_star[loc,tpt]
        }
    }
}


xDim = dim(X)
zDim = c(dim(Z_ar)[1], prod(dim(Z_ar)[2:3]))

S0 = S[,1]
E0 = E[,1]
I0 = I[,1]
R0 = R[,1]

S_star0 = rep(0, nrow(S))
E_star0 = rep(0, nrow(S))
I_star0 = rep(0, nrow(S))
R_star0 = rep(0, nrow(S))

#S = S[,2:ncol(S)]
#E = E[,2:ncol(S)]
#I = I[,2:ncol(S)]
#R = R[,2:ncol(S)]

#S_star = S_star[,2:ncol(S_star)]
#E_star = E_star[,2:ncol(S_star)]
#I_star = I_star[,2:ncol(S_star)]
#R_star = R_star[,2:ncol(S_star)]

compMatDim = c(nrow(S), ncol(S))


DM = as.numeric(neighborhood)

rho = 0.6

p_ei = 0.8
p_ir = 0.6
p_rs = rep(0.1, ncol(S))

beta = c(1,0.0,0.0)

outFileName = "./chainOutput.txt"
# beta, rho, p_se, p_ei, p_ir,p_rs,S*,E*,I*,R*
logFileList = c(1,1,0,1,1,1,0,0,0,0)
iterationStride = 1
res = spatialSEIRInit(compMatDim,xDim,
                      zDim,S0,
                      E0,I0,
                      R0,S_star0,
                      E_star0,I_star0,
                      R_star0,S_star,
                      E_star,I_star,
                      R_star,X,
                      Z,DM,
                      rho,beta,
                      p_ei,p_ir,
                      p_rs,N,outFileName, logFileList, iterationStride)


