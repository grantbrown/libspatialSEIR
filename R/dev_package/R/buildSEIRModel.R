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

buildSEIRModel = function(dataModel,
                          distanceModel,
                          exposureModel,
                          reinfectionModel,
                          samplingControl)
{

    interface = new( spatialSEIRModel(dataModel,
                                      distanceModel,
                                      reinfectionModel,
                                      samplingControl))    
}
