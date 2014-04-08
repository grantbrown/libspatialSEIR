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

#####################################
# Part 2: Format for libSpatialSEIR #
#####################################

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


