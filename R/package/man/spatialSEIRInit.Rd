\name{spatialSEIRInit}
\alias{spatialSEIRInit}

\title{
    Initialization function, passing data, instructions, and control to the C++ library. 
}
\description{
    Initialization function, passing data, instructions, and control to the C++ library. 
    Exported by the R library for testing and development purposes, a much friendlier 
    final library wrapper will be provided when the C++ library actually works. 
}
\usage{
    spatialSEIRInit(compMatDim,xDim,zDim,S0_,E0_,I0_,R0_,Sstar0,Estar0,Istar0,
                    Rstar0,Sstar, Estar,Istar,Rstar,X_,Z_,DistMat_,rho_,gamma_,
                    priorAlpha_gamma_, priorBeta_gamma_,priorAlpha_pEI_,priorBeta_pEI_,
                    priorAlpha_pIR_,priorBeta_pIR_,beta_,p_ei_,p_ir_,p_rs_,N_,outFile,
                    logVarList,iterationStride,verboseFlag,debugFlag,sliceWidths)
}
\arguments{
    \item{compMatDim}{
        Dimension (row column) of the compartmental model matrix object to be populated by Istar.
    }
    \item{xDim}{
        Dimension (row column) of the time fixed covariate matrix X
    }
    \item{zDim}{
        Dimension (row column) of the time varying covariate matrix Z
    }
    \item{S0_}{
        Initial vector for S
    }
    \item{E0_}{
        Initial vector for E
    }
    \item{I0_}{
        Initial vector for I
    }
    \item{R0_}{
        Initial vector for R
    }
    \item{Sstar0}{
        Initial R to S transition vector
    }
    \item{Estar0}{
        Initial S to E transition vector
    }
    \item{Istar0}{
        Initial E to I transition vector
    }
    \item{Rstar0}{
        Initial I to R transition vector
    }
    \item{Sstar}{
        Observed infection data
    }
    \item{Estar}{
        Observed infection data
    }
    \item{Istar}{
        Observed infection data
    }
    \item{Rstar}{
        Observed infection data
    }
    \item{X_}{
        Covariate data
    }
    \item{Z_}{
        Time varying covariate data
    }
    \item{DistMat_}{
        Raw distance matrix. 
    }
    \item{rho_}{
        Initial value for spatial dependence parameter.
    }
    \item{gamma_}{
        Vector of baseline infection parameters. 
    }
    \item{priorAlpha_gamma_}{
        Prior alpha value for gamma.
    }
    \item{priorBeta_gamma_}{
        Prior beta value for gamma.
    }
    \item{priorAlpha_pEI_}{
        Prior alpha value for pEI.
    }
    \item{priorBeta_pEI_}{
        Prior beta value for pEI.
    }

    \item{priorAlpha_pIR_}{
        Prior alpha value for pIR.
    }
    \item{priorBeta_pIR_}{
        Prior beta value for pIR.
    }
    \item{beta_}{
        Initial vector for exposure probability
    }
    \item{p_ei_}{
        Initial value for exposed to infectious probability.
    }
    \item{p_ir_}{
        Initial value for infectious to recovered probability.
    }
    \item{p_rs_}{
        Initial vector for recovered to susceptible probability.
    }
    \item{N_}{
        Population sizes at each spatial location (currently fixed)
    }
    \item{outFile}{
        Path to file where chain output will be stored. 
    }
    \item{logVarList}{
        An integer vector, which must currently be of length 10. Nonzero entries indicate the the corresponding set of parameters should be saved to the output file specified in outFile. The entries correspond to: beta, rho, p_se, p_ei, p_ir, p_rs, S*, E*, I*, and R*.   
    }
    \item{iterationStride}{
        Thinning to be used when saving iterations to disk. A stride of 1 will save all samples. 
    }
    \item{verboseFlag}{
        Show verbose output?
    }
    \item{debugFlag}{
        Show debug level output?
    }
    \item{sliceWidths}{
        Slice sampling tuning parameters for the slice sampling variables.
    }








}
\author{
    Grant Brown
}
\keyword{ utilities }

