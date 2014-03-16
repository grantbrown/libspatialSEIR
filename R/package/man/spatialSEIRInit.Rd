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
    spatialSEIRInit(compMatDim,xDim,zDim,S0_,E0_,I0_,R0_,Sstar0,Estar0,Istar0,Rstar0,Sstar, Estar,Istar,Rstar,X_,Z_)
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

}
\author{
    Grant Brown
}
\keyword{ utilities }

