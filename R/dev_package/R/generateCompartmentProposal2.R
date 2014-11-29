generateCompartmentProposal2 = function(data, N, S0=NA, E0=NA, I0=NA, reinfection=FALSE, dataType = c("I_star", "R_star"))
{
    
    data = as.matrix(data)
    N = as.matrix(N)
    if (!all(dim(N) == dim(data)))
    {
        stop("Error in generateCompartmentProposal2: N and data must be same dimension.")
    }
    if (length(S0) != length(E0) || length(E0) != length(I0) || length(I0) != ncol(data))
    {
        stop("Initial values imply different number of locations from data.")
    }
    if (any(c(is.na(S0), is.na(E0), is.na(I0))))
    {
        stop("None of: S0, E0, I0 may be missing.")
    }

    R0 = N[1,] - S0 - E0 - I0


    S = matrix(S0, nrow = nrow(data), ncol = length(S0), byrow=TRUE)
    E = matrix(E0, nrow = nrow(data), ncol = length(S0), byrow=TRUE)
    I = matrix(I0, nrow = nrow(data), ncol = length(S0), byrow=TRUE)
    R = matrix(R0, nrow = nrow(data), ncol = length(S0), byrow=TRUE)

    nTpt = nrow(S)
    nLoc = ncol(S)
    if (dataType[1] == "I_star")
    {
        E_star = rbind(data[2:nrow(data),,drop=FALSE],rep(0, ncol(S)))
        I_star = data
        R_star = rbind(rep(0, ncol(S)), data[1:(nrow(data)-1),,drop=FALSE]) 
        if (reinfection)
        {
            S_star = rbind(rep(0, ncol(S)), rep(0, ncol(S)),data[1:(nrow(data)-2),,drop=FALSE])
        }
        else
        {
            S_star = E_star*0
        }
    }
    else
    {
        I_star = rbind(data[2:nrow(data),],rep(0, ncol(S)))
        R_star = data
        if (reinfection)
        {
            S_star = rbind(rep(0, ncol(S)), rep(0, ncol(S)),data[1:(nrow(data)-2),,drop=FALSE])
        }
        else
        {
            S_star = I_star*0
        }
        E_star = rbind(rep(0, ncol(S)), data[1:(nrow(data)-1),,drop=FALSE]) 
    }

    S = S + apply(S_star, 2, cumsum) - apply(E_star, 2, cumsum)
    E = E + apply(E_star, 2, cumsum) - apply(I_star, 2, cumsum)
    I = I + apply(I_star, 2, cumsum) - apply(R_star, 2, cumsum)
    R = R + apply(R_star, 2, cumsum) - apply(S_star, 2, cumsum)

    if (any(S<0) && !reinfection)
    {
        cat("Error: invalid compartment generated. Is this data supposed to incorporate reinfection?\n")
        # Return anyway for debugging purposes
    }
    else if (any(S<0) || any(E <0) || any(I<0) || any(R<0))
    {
        cat("Error: invalid compartment generated.\n")
        # Return anyway for debugging purposes
    }

    return(list(S0=S0,
                E0=E0,
                I0=I0,
                R0=R0,
                S_star=S_star,
                E_star=E_star,
                I_star=I_star,
                R_star=R_star,
                S=S,
                E=E,
                I=I,
                R=R))
}



