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

gen_covariates = function(dcm, nweek = 52, nyear =5)
{
    cat("Generating covariates...\n")
    X = cbind(1,matrix(rnorm(3*(ncol(dcm))), ncol = 3), 
            rbinom(ncol(dcm), 1, 0.2)
            #, 
            #diag(rep(1, ncol(dcm))))
            )

    # Month, Temperature by year
    Z = array(0, dim = c(nrow(X),4,nweek, nyear))

    for (i in 1:nyear)
    {
        for (j in 1:nweek)
        {
            # Get actual temperature values here
            temp_mean = 2*sqrt(-(1:52-26.5)^2 + 800)
            Z[,,j,i] = cbind(i,(j-26.5),(j-26.5)^2, rnorm(nrow(X), temp_mean, 4)) 
        }
    }


    #true_fixed_beta = c(-1.2,-0.03,0,-0.025, -0.02, rnorm(ncol(dcm),0, 0.01))
    true_fixed_beta = c(-1.2,-0.03,0,-0.025, -0.02)

    true_time_varying_beta = c(0.004,0.0006, 0.001, -0.0004)
    
    true_fixed_eta = X %*% true_fixed_beta
    true_time_varying_eta = array(0, dim = c(nrow(X), nweek, nyear))

    for(i in 1:nyear)
    {
        for (j in 1:nweek)
        {
            true_time_varying_eta[,j,i] = Z[,,j,i] %*% true_time_varying_beta
        }
    }
    return(list("true_time_varying_eta"=true_time_varying_eta, 
                "true_fixed_eta"=true_fixed_eta, 
                "true_time_varying_beta" = true_time_varying_beta, 
                "true_fixed_beta" = true_fixed_beta, 
                "X" = X, 
                "Z" = Z))
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

main_sim = function(pop, dcm, X, Z, true_fixed_eta, true_time_varying_eta, true_rho = 0.8, nweek = 52, nyear = 5)
{
    cat("Running main simulation...\n")

    # TODO: expand to make N change by year
    N = pop[,2]

    S = array(0, dim = c(length(N), nweek, nyear))
    E = array(0, dim = c(length(N), nweek, nyear))
    I = array(0, dim = c(length(N), nweek, nyear))
    R = array(0, dim = c(length(N), nweek, nyear))

    S_star_output = array(0, dim = c(length(N), nweek, nyear))
    E_star_output = array(0, dim = c(length(N), nweek, nyear))
    I_star_output = array(0, dim = c(length(N), nweek, nyear))
    R_star_output = array(0, dim = c(length(N), nweek, nyear))

    S0 = floor(0.9*N)
    E0 = floor(0.005*N)
    I0 = rbinom(rep(1, length(N)), N, 0.05)
    R0 = N - S0 - E0 - I0
    R0 = ifelse(R0 >= 0, R0, 0)

    # Parameters:

    p_i = 0.8
    p_r = 0.3
    p_s = c(rep(0.005, 8), rep(0.05, 4), rep(0.1, 8), rep(0.4,4), rep(0.5,3),rep(0.925,3),rep(0.1,4),rep(0.005, 18))

    for (year in 1:nyear)
    {
        cat(paste(year, "\n", sep = ""))
        for (week in 1:nweek)
        {
            cat(".")
            if (week == 1 & year == 1)
            {
                # Initial Case
                # Compute infection probability

                mu = exp(true_fixed_eta + true_time_varying_eta[,week,year])

                # Calculate Infectious Ratio

                delta = I0/N
                dmatrix = (1/sqrt(dcm))*true_rho
                diag(dmatrix) = 0
                p = 1-exp(-delta*mu - as.numeric(dmatrix %*% (delta*mu)))

                if (any(is.na(p)))
                {
                    stop("Invalid P")
                }

                E_star0 = rbinom(rep(1, length(S0)),S0, p) 
                I_star0 = rbinom(rep(1, length(E0)),E0, p_i)
                R_star0 = rbinom(rep(1, length(I0)),I0, p_r)
                S_star0 = rbinom(rep(1, length(R0)),R0, p_s[week])

                if (any(is.na(E_star0)))
                {
                    stop("Invalid I star 0")
                }
                if (any(is.na(I_star0)))
                {
                    stop("Invalid I star 0")
                }
                if (any(is.na(R_star0)))
                {
                    stop("Invalid R star 0")
                }
                if (any(is.na(S_star0)))
                {
                    stop("Invalid S star 0")
                }

                S[,week,year] = S0 - E_star0 + S_star0
                E[,week,year] = E0 + E_star0 - I_star0
                I[,week,year] = I0 + I_star0 - R_star0
                R[,week,year] = R0 + R_star0 - S_star0
            }
            else
            {
                # Standard Case          
                # Compute infection probability

                if (week == 1){mr_week = 52; mr_year = year-1}
                else {mr_week = week-1; mr_year = year} 

                mu = exp(true_fixed_eta + true_time_varying_eta[,mr_week,mr_year])
                # Calculate Infectious Ratio
                delta = I[,mr_week,mr_year]/N
                dmatrix = delta*(1/dcm)*true_rho
                diag(dmatrix) = 0
                p = 1-exp(-delta*mu - as.numeric(dmatrix %*% (delta*mu)))

                if (any(is.na(p)))
                {
                    stop("Invalid P")
                }

                E_star = rbinom(rep(1, length(S0)),S[,mr_week, mr_year], p)
                I_star = rbinom(rep(1, length(E0)),E[,mr_week,mr_year], p_i)
                R_star = rbinom(rep(1, length(I0)),I[,mr_week,mr_year], p_r)
                S_star = rbinom(rep(1, length(R0)),R[,mr_week,mr_year], p_s[week])
                
                S_star_output[,mr_week,mr_year] = S_star
                E_star_output[,mr_week,mr_year] = E_star
                I_star_output[,mr_week,mr_year] = I_star
                R_star_output[,mr_week,mr_year] = R_star

                S[,week,year] = S[,mr_week, mr_year] + S_star - E_star
                E[,week,year] = E[,mr_week,mr_year] + E_star - I_star
                I[,week,year] = I[,mr_week,mr_year] + I_star - R_star
                R[,week,year] = R[,mr_week,mr_year] + R_star - S_star
            }

        }
    }
    return(list("S" = S, "E" = E, "I"=I, "R"=R, 
                "S_star" = S_star_output,
                "E_star" = E_star_output,
                "I_star" = I_star_output,
                "R_star" = R_star_output,
                "S0"=S0,
                "E0"=E0,
                "I0"=I0,
                "R0"=R0,
                "S_star0"=S_star0,
                "E_star0"=E_star0,
                "I_star0"=I_star0,
                "R_star0"=R_star0
                ))
}

plot_yearly_trends = function(SEIR_Object, main = "SEIR Plot")
{
    SEIR_dim = dim(SEIR_Object)
    plot(0,0,xlim = c(0,52*SEIR_dim[3]), ylim = c(0,max(SEIR_Object)), main = main, type = "n")

    yrs = SEIR_Object[,,1]
    for (year in 1:SEIR_dim[3])
    {
        yrs = cbind(yrs, SEIR_Object[,,year])
    }


    for (spat_loc in 1:SEIR_dim[1])
    {
        lines(x = seq(1,ncol(yrs)), y = yrs[spat_loc,], #col = col2rgb("grey", alpha = 0.5))
        col = "black")
    }
    #lines(x = matrix(seq(1, ncol(yrs)), nrow = nrow(yrs), ncol = ncol(yrs), byrow = TRUE), y = yrs)
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
    # Population Data
    pop = data_list$pop[keepLocations,] 

    covariates = gen_covariates(dcm)
    X = covariates$X
    Z = covariates$Z

    true_fixed_eta = covariates$true_fixed_eta
    true_time_varying_eta = covariates$true_time_varying_eta
    
    true_fixed_beta = covariates$true_fixed_beta
    true_time_varying_beta = covariates$true_time_varying_beta

    sim_results = main_sim(pop,dcm,X,Z,true_fixed_eta, true_time_varying_eta)

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

    data_list = list("mcl"=mcl, "dcm"=dcm, "pop"=pop)
    covariates = list("X"= X, "Z"=Z, "true_fixed_eta"=true_fixed_eta, "true_time_varying_eta"=true_time_varying_eta,"true_fixed_beta" = true_fixed_beta, "true_time_varying_beta" = true_time_varying_beta)
    save(data_list, covariates, sim_results, file = "./SimulationObjects.robj")
    cat("Done.\n")
}
