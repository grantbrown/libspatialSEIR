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

#I_star = Y[,2:ncol(Y)]
#I_star0 = Y[,1] 

#R_star = rbinom(Y, 0.8)[,2:ncol(Y)]
#R_star0 = rbinom(Y[,1], 0.8)

#E_star = Y[,1:(ncol(Y)-1)]
#E_star0 = Y[,ncol(Y)]  

E_star = Y[,c(2:(ncol(Y)-1),1)]
I_star = Y

#R_star = matrix(rbinom(rep(1, prod(dim(Y)),Y, 0.9)), nrow = nrow(Y), ncol = ncol(Y))
#E_star = cbind(Y[,2:(ncol(Y))], Y[,1])
#I = t(apply(I_star, 1, cumsum) - apply(R_star,1,cumsum))
#E = (Y[,1] + t((apply(E_star,1,cumsum))) - t(apply(I_star,1,cumsum)))

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
