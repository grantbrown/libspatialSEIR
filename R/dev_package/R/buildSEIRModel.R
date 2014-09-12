buildDistanceModel = function(distanceList, 
                              scaleMode = c("none","rowscale","invsqrt"))
{
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

buildDataModel = function(Y, modelType = c("identity", "overdispersion"))
{
    return(new(dataModel(Y, modelType)))
}

buildReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, priorPrecision = NA)
{
    integerMode = ifelse(reinfectMode[1] == "SEIR", 3, 
                  ifelse(reinfectMode[1] == "SEIRS", 1, 
                  ifelse(reinfectMode[1] == "Fixed", 2, NA)))
    if (is.na(integerMode))
    {
        stop(paste("Invalid mode: ", reinfectMode[1], sep = ""))
    }
    if (integerMode != 3 && (is.na(X_prs)))
    {
        stop("If reinfection mode is not SEIR, X_prs must be supplied.")
    }
    if (integerMode == 1 && (is.na(priorPrecision)))
    {
        stop("If reinfection parameters are going to be estimated, priorPrecision must be specified.")
    }

    reinfectionmod = new(reinfectionModel(integerMode));
    if (integerMode != 3)
    {
        if (all(is.na(priorPrecision)))
        {
           priorPrecision = 0.1     
        }
        reinfectionmod$buildReinfectionModel(X_prs, priorPrecision);
    }
    reinfectionmod
}

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

buildSEIRModel = function(dataModelInstance,
                          exposureModelInstance,
                          reinfectionModelInstance,
                          distanceModelInstance,
                          transitionPriorsInstance,
                          initialValueContainer,
                          samplingControlInstance)
{

    interface = new( spatialSEIRModel )
    interace$buildSpatialSEIRModel(dataModelInstance,
                                   exposureModelInstance,
                                   reinfectionModelInstance,
                                   distanceModelInstance,
                                   transitionPriorsInstance,
                                   initialValueContainer,
                                   samplingControlInstance)
    interface
}
