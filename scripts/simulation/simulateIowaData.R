read_data = function()
{
    cat("Reading population data...\n")
    PopDat = read.table("../scriptData/census/ProcessedPopulationData.csv",
                        sep = ";", head = TRUE)
    cat("Reading distance data...\n")
    DistanceDat = read.table("../scriptData/OSM/DistancesProcessed.delim",
                            sep = ";", head = FALSE)
    colnames(DistanceDat) = c("City1", "City2", "Secs", "Mi")
    DistanceDat[,1][DistanceDat[,1] == "De Witt"] = "DeWitt"
    DistanceDat[,2][DistanceDat[,2] == "De Witt"] = "DeWitt"
    PopDat[,1] = as.character(PopDat[,1])
    PopDat[,1][PopDat[,1] == "De Witt"] = "DeWitt"
    PopDat = PopDat[order(PopDat[,1]),]

    DistanceDat = as.data.frame(DistanceDat)

    list("pop"=PopDat, "dist"=DistanceDat) 
}

gen_dist_mat = function(DistanceDat, PopDat)
{

    if (file.exists("./distmat.robj"))
    {
        load("./distmat.robj")
        return(distmatlist)
    }

    u1= unique(c(as.character(DistanceDat$City1), as.character(DistanceDat$City2)))
    u2  = unique(as.character(PopDat$City))
    masterCityList = intersect(u1, u2)
    masterCityList = masterCityList[order(masterCityList)]

    indexDistance1 = function(drow)
    {
        return(which(masterCityList == drow[1]))
    }    
    indexDistance2 = function(drow)
    {
        return(which(masterCityList == drow[2]))
    }
    cat("Indexing distance measures...\n")
    idx1 = as.numeric(apply(DistanceDat, 1, indexDistance1))
    idx2 = as.numeric(apply(DistanceDat, 1, indexDistance2))
    # Remove cities not in master city list:
    keep_indices = !(is.na(idx1) | is.na(idx2))
    idx1 = idx1[keep_indices]
    idx2 = idx2[keep_indices]
 
    dcm = diag(rep(0, length(masterCityList)))
    cat("Populating distance matrix...\n")
    dcm[cbind(idx1, idx2)] = matrix(DistanceDat$Secs[keep_indices])
    # Fix triangularity
    tdcm = t(dcm)
    gtIdx = tdcm > dcm
    dcm[gtIdx] = tdcm[gtIdx]

    colnames(dcm) = masterCityList
    rownames(dcm) = masterCityList

    distmatlist = list("mcl"=masterCityList,"idx1"=idx1, "idx2"=idx2, "dcm"=dcm, "keep" = keep_indices)
    save(distmatlist, file = "distmat.robj")
    distmatlist
}

protected_log = function(x, min_x = 0.001)
{
    out = log(x)
    return(ifelse(out == -Inf, log(min_x), out))
}

logit = function(x)
{
    return(log(x/(1-x)))
}
ilogit = function(x)
{
    ex = exp(x)
    return(ex/(1+ex))
}





main_sim = function(dcm, pop, nTptPerYear = 12, nyear =5)
{
    cat("Generating covariates...\n")

    maxTpt = nTptPerYear*nyear
    # For simplicity just use intercept for now. 
    X = matrix(1, ncol = 1, nrow = nrow(dcm))
    Z = cbind(seq(1,maxTpt), model.matrix(~as.factor(rep(1:12, nyear)))[,2:nTptPerYear])
    Z = Z[rep(1:nrow(Z), each = nrow(dcm)),] # Same time varying covariates for all spatial locations 

    X_prs = cbind(1, 
                  sin((1:maxTpt)/nTptPerYear*2*pi),
                  cos((1:maxTpt)/nTptPerYear*2*pi))

    trueBetaSEFixed = c(-0.2)
    trueBetaSEVarying = c(0.0001, 0.02, 0.03, 0.05, 0.06, 
                          0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1)*5

    trueGamma = rep(0.0, maxTpt)
    trueRho = 0.05
    rho = trueRho

    true_fixed_eta = X %*% trueBetaSEFixed
    true_time_varying_eta = Z %*% trueBetaSEVarying

    eta_se = matrix((as.numeric(true_fixed_eta) + 
                     as.numeric(true_time_varying_eta)), 
                     nrow = maxTpt, ncol = ncol(dcm), byrow = TRUE)

    p_se = matrix(0.0, nrow = maxTpt, ncol = ncol(dcm)) 
    p_ei = 0.9
    p_ir = 0.9
    trueBetaRS = c(2.5, -1, 0.5)
    eta_rs = X_prs %*% trueBetaRS    
    p_rs = exp(-eta_rs)

    cat("Running main simulation...\n")

    # TODO: expand to make N change by year
    N = pop[,2]
    S = matrix(0,nrow = maxTpt, ncol = length(N))
    E = S 
    I = S
    R = S

    S_star = S
    E_star = S
    I_star = S 
    R_star = S

    E0 = 0*N
    I0 = floor(0.001*N)
    R0 = floor(0.001*N) 
    S0 = N - R0 - E0 - I0

    S_star0 = 0*N
    E_star0 = I0
    I_star0 = 0*N
    R_star0 = floor(0.99*I0)


    # Parameters:



    for (i in 1:maxTpt)
    {
        if (i == 1)
        {
            S[i,] = S0 + S_star0 - E_star0
            E[i,] = E0 + E_star0 - I_star0
            I[i,] = I0 + I_star0 - R_star0
            R[i,] = R0 + R_star0 - S_star0
            etaVal = (I[i,]/N)*exp(eta_se[i,])
            etaVal = etaVal + rho*(dcm %*% etaVal) + trueGamma[i] 
            p_se[i,] = 1-exp(-etaVal)

            S_star[i,] = rbinom(rep(1,ncol(S)), R[i,], p_rs[i])
            E_star[i,] = rbinom(rep(1,ncol(S)), S[i,], p_se[i,])
            I_star[i,] = rbinom(rep(1,ncol(S)), E[i,], p_ei)
            R_star[i,] = rbinom(rep(1,ncol(S)), I[i,], p_ir)
        }
        else
        {
            S[i,] = S[i-1,] + S_star[i-1,] - E_star[i-1,]
            E[i,] = E[i-1,] + E_star[i-1,] - I_star[i-1,]
            I[i,] = I[i-1,] + I_star[i-1,] - R_star[i-1,]
            R[i,] = R[i-1,] + R_star[i-1,] - S_star[i-1,]

            etaVal = (I[i,]/N)*exp(eta_se[i,])
            etaVal = etaVal + rho*(dcm %*% etaVal) + trueGamma[i] 
            p_se[i,] = 1-exp(-etaVal)

            S_star[i,] = rbinom(rep(1,ncol(S)), R[i,], p_rs[i])
            E_star[i,] = rbinom(rep(1,ncol(S)), S[i,], p_se[i,])
            I_star[i,] = rbinom(rep(1,ncol(S)), E[i,], p_ei)
            R_star[i,] = rbinom(rep(1,ncol(S)), I[i,], p_ir)

        }
    }

    return(list("S" = S, "E" = E, "I"=I, "R"=R, 
                "S_star" = S_star,
               "E_star" = E_star,
                "I_star" = I_star,
                "R_star" = R_star,
                "S0"=S0,
                "E0"=E0,
                "I0"=I0,
                "R0"=R0,
                "S_star0"=S_star0,
                "E_star0"=E_star0,
                "I_star0"=I_star0,
                "R_star0"=R_star0,
                "gamma" = trueGamma,
                "true_fixed_beta" = trueBetaSEFixed,
                "true_time_varying_beta" = trueBetaSEVarying,
                "true_beta_pRS"=trueBetaRS,
                "X" = X,
                "Z" = Z,
                "X_prs" = X_prs,
                "p_ei" = p_ei,
                "p_ir" = p_ir,
                "p_se" = p_se
                ))
}

control_code = function(Nlocations = 20, plots = FALSE)
{
    data_list = read_data()
    distance_stuff = gen_dist_mat(data_list[["dist"]], data_list[["pop"]])

    keepLocations = sample(1:ncol(distance_stuff[["dcm"]]))[1:Nlocations]

    # Create Final Data Subsets
    # Master City List
    mcl=distance_stuff$mcl[keepLocations]
    # Distance Matrix
    dcm=distance_stuff$dcm[keepLocations,keepLocations] 
    sDcm = 1/sqrt(dcm)
    sDcm = ifelse(dcm < 4*60*60, sDcm, 0)
    dcm = ifelse(dcm < 4*60*60, dcm, 0)
    diag(sDcm) = 0
    # Population Data
    pop = data_list$pop[keepLocations,] 

    sim_results = main_sim(sDcm,pop)

    if (plots)
    {
        par(mfrow = c(4,2))
        plot_yearly_trends(sim_results$S, main ="S")
        plot_yearly_trends(sim_results$S_star, main = "S*")
        plot_yearly_trends(sim_results$E, main = "E")
        plot_yearly_trends(sim_results$E_star, main = "E*")
        plot_yearly_trends(sim_results$I, main = "I")
        plot_yearly_trends(sim_results$I_star, main = "I*")
        plot_yearly_trends(sim_results$R, main = "R")
        plot_yearly_trends(sim_results$R_star, main = "R*")
        cat("\n Diagnostics:\n ")
        cat("S_star0:\n", sim_results$S_star0, "\n")
        cat("E_star0:\n", sim_results$E_star0, "\n")
        cat("I_star0:\n", sim_results$I_star0,"\n")
        cat("R_star0:\n", sim_results$R_star0,"\n")
    }

    data_list = list("mcl"=mcl, "dcm"=dcm, "pop"=pop, "sDcm"=sDcm) 
    save(data_list, sim_results, file = "./SimulationObjects.robj")
    cat("Done.\n")
}






