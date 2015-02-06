## distanceModel module helper function
buildDistanceModel = function(distanceList, 
                              scaleMode = c("none","rowscale","invsqrt"),
                              priorAlpha=1.0,
                              priorBeta=1.0)
{
    scaleMode = scaleMode[1]
    rowScale = function(mat)
    {
        mat/matrix(apply(mat,1,sum), nrow = nrow(mat), ncol = ncol(mat))
    }

    invSqrt = function(mat)
    {
        matrix(ifelse(mat == 0, 0, 1/sqrt(mat)), nrow = nrow(mat), ncol = ncol(mat))
    }


    # Check if we've got a single matrix
    if (class(distanceList) == "matrix")
    {
        distanceList = list(distanceList)
    }
    # Make sure that at this point we have a list of matrices. 

    if (class(distanceList) != "list")
    {
        stop("Error: distanceList must be a list of matrices.")
    }

    # Check for valid matrices
    distanceDim = NA
    for (i in 1:length(distanceList))
    {
        if (class(distanceList[[i]]) != "matrix")
        {
            stop("Distance metrics must be matrices.")
        }
        distanceDim = dim(distanceList[[1]])
        newDim = dim(distanceList[[i]])
        if (any(newDim != distanceDim) || (newDim[1] != newDim[2]))
        {
            stop("Distance matrices must be square and of the same dimension.")
        }
        if (any(is.na(distanceList[[i]])))
        {
            stop("NA's are not allowed in distance matrices.")
        }
        if (scaleMode == "rowscale")
        {
            distanceList[[i]] = rowScale(distanceList[[i]])
        }
        if (scaleMode == "invsqrt")
        {
            distanceList[[i]] = invSqrt(distanceList[[i]])
        }
    } 

    finalDistanceModel = new( distanceModel )
    for (i in 1:length(distanceList))
    {
        finalDistanceModel$addDistanceMatrix(distanceList[[i]])
    }
    finalDistanceModel$setPriorParameters(priorAlpha, priorBeta);
    finalDistanceModel
}

# dataModel module helper function
buildDataModel = function(Y, type = c("identity", "overdispersion"), compartment = c("I_star", "R_star"),params=NA)
{
    if (class(Y) != "matrix")
    {
        Y = as.matrix(Y)
    }
    type = type[1]
    if (length(params) == 1 && is.na(params) && type != "identity")
    {
        stop("Non-identity data model selected without specifying parameters.")
    }
    else if (length(params) ==1 && is.na(params))
    {
        return(new(dataModel, Y, type, compartment))
    }
    else if (class(params) != "numeric" || length(params) != 2)
    {
        stop("Non identy data model currently requires the params argument to be a numeric vector of length 2") 
    }
    outModel = new(dataModel, Y, type, compartment)
    outModel$setOverdispersionParameters(params[1], params[2], rgamma(1, params[1], params[2]))
    outModel
}

# reinfectionModel module helper function
buildReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, 
                                 betaPrs=NA , priorPrecision = NA, priorMean = NA)
{
    integerMode = ifelse(reinfectMode[1] == "SEIR", 3, 
                  ifelse(reinfectMode[1] == "SEIRS", 1, 
                  ifelse(reinfectMode[1] == "Fixed", 2, NA)))
    if (is.na(integerMode))
    {
        stop(paste("Invalid mode: ", reinfectMode[1], sep = ""))
    }
    if (integerMode != 3 && (is.na(X_prs) || is.na(betaPrs)))
    {
        stop("If reinfection mode is not SEIR, X_prs and betaPrs must be supplied.")
    }
    if (integerMode == 1 && (is.na(priorPrecision) || is.na(priorMean)))
    {
        stop("If reinfection parameters are going to be estimated, priorPrecision and priorMean must be specified.")
    }

    reinfectionmod = new(reinfectionModel, integerMode);
    if (integerMode != 3)
    {
        if (class(X_prs) != "matrix")
        {
            X_prs = as.matrix(X_prs)
        }
        if (all(is.na(priorPrecision)))
        {
           priorPrecision = rep(0.1, ncol(X_prs))     
        }
        else if (length(priorPrecision) == 1)
        {
            priorPrecision = rep(priorPrecision, ncol(X_prs))               
        }
        if (all(is.na(priorMean)))
        {
            priorMean = rep(0, ncol(X_prs))
        }
        else if (length(priorMean) == 1)
        {
            priorMean = rep(priorMean, ncol(X_prs))
        }
        reinfectionmod$buildReinfectionModel(X_prs, betaPrs, priorMean, priorPrecision);
    }
    reinfectionmod
}

# samplingControl module helper function
buildSamplingControl = function(verbose = FALSE, debug = FALSE, iterationStride = 100, steadyStateConstraintPrecision = -1, sliceWidths = rep(0.1, 11))
{
    samplingControlInstance = new ( samplingControl )
    samplingControlInstance$verbose = verbose
    samplingControlInstance$debug = debug
    samplingControlInstance$iterationStride = iterationStride
    samplingControlInstance$steadyStateConstraintPrecision = steadyStateConstraintPrecision
    samplingControlInstance$sliceWidths = sliceWidths
    samplingControlInstance        
}

# transitionPriors module helper functions
buildTransitionPriorsFromProbabilities = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    tp = new(transitionPriors)
    tp$setPriorsFromProbabilities(p_ei, p_ir, p_ei_ess, p_ir_ess)
    tp
}
buildTransitionPriorsManually = function(priorAlpha_gammaEI, priorBeta_gammaEI,
                                         priorAlpha_gammaIR, priorBeta_gammaIR)
{
     tp = new(transitionPriors)
     tp$setPriorsManually(priorAlpha_gammaEI, priorBeta_gammaEI,
                          priorAlpha_gammaIR, priorBeta_gammaIR)
     tp
}
buildUniformTransitionPriors = function()
{
    tp = new(transitionPriors)
    tp$setUniformPriors()
    tp
}

# exposureModel module helper function
buildExposureModel = function(X,nTpt, nLoc, beta=NA,betaPriorPrecision=NA,
                              betaPriorMean=NA,offset=NA)
{
    nBeta = ncol(X)
    if (class(X) != "matrix")
    {
        print("Warning: X should be a matrix.")
    }
    if (nrow(X) != nTpt * nLoc){
        stop("Invalid data dimensions: X should be a matrix composed of nLoc row-wise blocks of dimension nTpt*p.") 
    }
    if (length(beta) == 1 && is.na(beta))
    {
        print("Generating starting values for exposure parameters: may be unreasonable.")
        beta = rnorm(nBeta)
    }
    if (length(betaPriorPrecision) == 1 && is.na(betaPriorPrecision))
    {
        print("No prior precision specified, using zero.")
        betaPriorPrecision = rep(0.1, nBeta)
    }
    else if (length(betaPriorPrecision) == 1)
    {
        betaPriorPrecision = rep(betaPriorPrecision, nBeta)
    }
    if (length(betaPriorMean) == 1 && is.na(betaPriorMean))
    {
        print("No prior mean specified, using zero.")
        betaPriorMean = rep(0, nBeta)
    }
    else if (length(betaPriorMean) == 1)
    {
        betaPriorMean = rep(betaPriorMean, nBeta)
    }
    if (length(offset) == 1 && is.na(offset))
    {
        print("Assuming equally spaced count data.")
    }
    else if (length(offset) == 1){
        offset = rep(offset, nTpt)
    }
    else if (length(offset) != nTpt){
        stop("Offset must be the of length nTpt")
    }
    ExposureModel = new(exposureModel,X,nTpt,nLoc,beta,betaPriorMean,betaPriorPrecision)
    if (all(!is.na(offset)))
    {
        ExposureModel$offsets = offset
    }
    ExposureModel
}

# depricated exposure interace employing separate X and Z matrices. 
# will be removed in a future version of the spatialSEIR R package
buildExposureModel_depricated = function(X,Z=NA,beta=NA,betaPriorPrecision=NA,
                                         betaPriorMean=NA,offset=NA,nTpt=NA)
{
  hasZ = !(length(Z) == 1 && is.na(Z))
  nBeta = ncol(X) + ifelse(hasZ, ncol(Z), 0)
  if (class(X) != "matrix")
  {
    print("Warning: X should be a matrix.")
  }
  if (length(beta) == 1 && is.na(beta))
  {
    print("Generating starting values for exposure parameters: may be unreasonable.")
    beta = rnorm(nBeta)
  }
  if (length(betaPriorPrecision) == 1 && is.na(betaPriorPrecision))
  {
    print("No prior precision specified, using zero.")
    betaPriorPrecision = rep(0.1, nBeta)
  }
  else if (length(betaPriorPrecision) == 1)
  {
    betaPriorPrecision = rep(betaPriorPrecision, nBeta)
  }
  if (length(betaPriorMean) == 1 && is.na(betaPriorMean))
  {
    print("No prior mean specified, using zero.")
    betaPriorMean = rep(0, nBeta)
  }
  else if (length(betaPriorMean) == 1)
  {
    betaPriorMean = rep(betaPriorMean, nBeta)
  }
  if (length(offset) == 1 && is.na(offset))
  {
    print("Assuming equally spaced count data.")
  }
  if (length(Z) == 1 && is.na(Z) && length(nTpt) == 1 && is.na(nTpt))
  {
    stop("If time varying covariate matrix Z is not provided, nTpt must be specified.")
  }
  if (length(Z) == 1 && is.na(Z))
  {
    nLoc = nrow(X)
    X = X[rep(1:nrow(X), each = nTpt),]
    ExposureModel = new(exposureModel,X,nTpt,nLoc,beta,betaPriorMean,betaPriorPrecision) 
  }
  else
  {
    nLoc = nrow(X)
    nTpt = floor(nrow(Z)/nrow(X))
    X = X[rep(1:nrow(X), each = nTpt),]
    X = cbind(X, Z)
    ExposureModel = new(exposureModel,X,nTpt,nLoc,beta,betaPriorMean,betaPriorPrecision)
  }
  if (all(!is.na(offset)))
  {
    ExposureModel$offsets = offset
  }
  ExposureModel
}

# initialValueContainer module helper function
buildInitialValueContainer = function(data, N, S0=NA, E0=NA, I0=NA, reinfection=FALSE, dataType=c("I_star", "R_star"))
{
    # get rid of any data frame nonsense
    S0 = as.numeric(S0);
    E0 = as.numeric(E0);
    I0 = as.numeric(I0);
    if (class(data) != "matrix")
    {
        data = as.matrix(data)
    }
    if (class(N) != "matrix")
    {
        N = as.matrix(N)
    }

    proposal = generateCompartmentProposal2(data, N, S0, E0, I0, reinfection, dataType)
    InitialValueContainer = new(initialValueContainer)
    InitialValueContainer$setInitialValues(proposal$S0, proposal$E0, proposal$I0, proposal$R0,
                                           proposal$S_star, proposal$E_star, proposal$I_star, proposal$R_star,
                                           N)
    InitialValueContainer
}

buildInitialValueContainerManually = function(S0=NA, E0=NA, I0=NA, R0=NA, S_star=NA, E_star=NA, I_star=NA, R_star=NA, N=NA)
{
    if (any(is.na(S0)) || any(is.na(E0)) || any(is.na(I0)) || any(is.na(R0)) 
        || any(is.na(S_star)) || any(is.na(E_star)) || any(is.na(I_star)) 
        || any(is.na(R_star)) || any(is.na(N)))
    {
        stop("Must specify all of S0, E0, I0, R0, S_star, E_star, I_star, R_star, N")
    }
    InitialValueContainer = new(initialValueContainer)
    InitialValueContainer$setInitialValues(S0, E0, I0, R0, S_star, E_star, I_star, R_star, N)
    InitialValueContainer
}


# SEIRModel module helper function
buildSEIRModel = function(outFileName,
                          dataModelInstance,
                          exposureModelInstance,
                          reinfectionModelInstance,
                          distanceModelInstance,
                          transitionPriorsInstance,
                          initialValueContainer,
                          samplingControlInstance)
{

    interface = new( spatialSEIRModel,outFileName )
    interface$buildSpatialSEIRModel(dataModelInstance,
                                   exposureModelInstance,
                                   reinfectionModelInstance,
                                   distanceModelInstance,
                                   transitionPriorsInstance,
                                   initialValueContainer,
                                   samplingControlInstance)
    interface
}
