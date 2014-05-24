generateCompartmentProposal = function(I_star, N, S0 = NA, E0 = NA, I0 = NA)
{
    if (dim(I_star)[1] <= 1)
    {
        stop("Invalid Compartment Dimension: Should be TxP")
    }
    singleLocation = ifelse(dim(I_star)[2] == 1, 1, 0) 
    if (any(dim(N) != dim(I_star)))
    {
        if ((dim(N)[0] == 2) && (dim(N)[2] == dim(I_star)[2]))
        {
            N = matrix(N, nrow = nrow(I), ncol = ncol(I), byrow = TRUE)
        }
        else
        {
            stop("N dimension does not match I_star dimension.")
        }
    }

    if (any(is.na(S0)))
    {
        S0 = floor(0.95*N[1,])
    }
    if (any(is.na(E0)))
    {
        E0 = I_star[1,]
    }
    if (any(is.na(I0)))
    {
        I0 = floor(I_star[1,]/2)
    }
   
    R0 = N[1,] - floor(S0 + E0 + I0) 
    S=E=I=R=S_star=E_star=R_star=I_star*0 

    for (tpt in 1:nrow(R))
    {
        if (tpt == 1)
        {
            S[1,] = S0
            E[1,] = E0
            I[1,] = I0
            R[1,] = R0

            S_star[1,] = rbinom(rep(1, length(S0)), R[1,], 0.05)
            E_star[1,] = rbinom(rep(1, length(S0)), I_star[2,], 1)
            R_star[1,] = rbinom(rep(1, length(S0)), I[1,], 0.8)
        }
        else
        {
            S[tpt,] = S[tpt-1,] + S_star[tpt-1,] - E_star[tpt-1]
            E[tpt,] = E[tpt-1,] + E_star[tpt-1,] - I_star[tpt-1]
            I[tpt,] = I[tpt-1,] + I_star[tpt-1,] - R_star[tpt-1]
            R[tpt,] = R[tpt-1,] + R_star[tpt-1,] - S_star[tpt-1]

            S_star[tpt,] = rbinom(rep(1, length(S0)), R[tpt,], 0.05)
            if (tpt != nrow(R))
            {
                E_star[tpt,] = I_star[tpt+1,]
            }
            else
            {
                E_star[tpt,] = E_star[tpt-1,]
            }
            R_star[tpt,] = rbinom(rep(1,length(S0)), I[tpt,],0.8)
        }
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



