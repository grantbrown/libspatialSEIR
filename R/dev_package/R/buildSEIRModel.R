## distanceModel module helper function
buildDistanceModel = function(distanceList, 
                              scaleMode = c("none","rowscale","invsqrt"))
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
    finalDistanceModel
}

# dataModel module helper function
buildDataModel = function(Y, type = c("identity", "overdispersion"),params=NA)
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
        return(new(dataModel, Y, type))
    }
    else if (class(params) != "numeric" || length(params) != 2)
    {
        stop("Non identy data model currently requires the params argument to be a numeric vector of length 2") 
    }
    outModel = new(dataModel, Y, type)
    outModel$setOverdispersionParameters(params[1], params[2], rgamma(1, params[1], params[2]))
    outModel
}

# reinfectionModel module helper function
buildReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, 
                                 betaPrs=NA , priorPrecision = NA)
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
    if (integerMode == 1 && (is.na(priorPrecision)))
    {
        stop("If reinfection parameters are going to be estimated, priorPrecision must be specified.")
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
           priorPrecision = 0.1     
        }
        reinfectionmod$buildReinfectionModel(X_prs, betaPrs, priorPrecision);
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
buildExposureModel = function(X,Z,beta,betaPriorPrecision,offset=NA)
{
    ExposureModel = new(exposureModel,X,Z,beta,betaPriorPrecision)
    if (!is.na(offset))
    {
        ExposureModel$offsets = offset
    }
    ExposureModel
}

# initialValueContainer module helper function

buildInitialValueContainer = function(I_star, N, S0=NA, E0=NA, I0=NA, reinfection =TRUE, p_ir = 0.9, p_rs = 0.05,ensureConstantInfectious = FALSE)
{
    if (class(I_star) != "matrix")
    {
        I_star = as.matrix(I_star)
    }
    if (class(N) != "matrix")
    {
        N = as.matrix(N)
    }


    proposal = generateCompartmentProposal(I_star, N, S0, E0, I0, reinfection, p_ir, p_rs,ensureConstantInfectious)
    InitialValueContainer = new(initialValueContainer)
    InitialValueContainer$setInitialValues(proposal$S0, proposal$E0, proposal$I0, proposal$R0,
                                           proposal$S_star, proposal$E_star, proposal$I_star, proposal$R_star,
                                           N)
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
