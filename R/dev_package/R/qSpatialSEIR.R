qSpatialSEIR <- function(formula, N, spatial.factor, distance.list=NA, verbose=TRUE, p_ei=NA, p_ir=NA, transition_ess=NA, chainFileName=NA, seed=NA, data, offset){
  # The following is based on the built in 'lm' function, 
  # and handles the formula and data provided by the user.
  call <- match.call()
  if (verbose){
    cat("Processing user input.\n")
    cat("Call: ")
    print(call)
    cat("\n")
  }
  if (any(is.na(p_ei)) | any(is.na(p_ir))){
    warning("Not specifying meaningful prior values for the E-I or I-R transition probabilities can result in model fitting problems.")
    if(any(is.na(p_ei))){
      p_ei = 0.5
      p_ei_ess = 1
    }
    if(any(is.na(p_ir))){
      p_ir = 0.5
      p_ir_ess = 1
    }
  }
  else{
    if (all(is.na(transition_ess))){
      p_ei_ess = 1000
      p_ir_ess = 1000
    }
    else{
      p_ei_ess = transition_ess 
      p_ir_ess = transition_ess 
    }
  }
  if (all(is.na(seed))){
    seed = as.numeric(Sys.time())
    set.seed(seed)
  }
  if (all(is.na(chainFileName))){
    chainFileName = paste("chainOutput_", seed, ".csv", sep = "")
  }
  if (missing(data))
  {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)  
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  cat("\n")
  I_star = model.response(mf, "numeric")
  mt <- attr(mf, "terms")
   
  if (is.empty.model(mt)){
    stop("Error, model formula was empty.")
  }  
  x <- model.matrix(mt, mf)

  # Process distance matrices and spatial.factor
  if (class(spatial.factor) != "factor"){
      spatial.factor = as.factor(spatial.factor)
  }

  if (all(is.na(distance.list)) || class(distance.list) != "list")
  {
      warning("Missing or invalid list of distance matrices. Assuming a single global spatial dependence term.")
      DM = (1-diag(length(levels(spatial.factor))))/(length(levels(spatial.factor)))
      distance.list = list(DM)
  }
  n.spatial.units = length(table(spatial.factor))
  for (dm in distance.list){
     if (class(dm) != "matrix" || nrow(dm) != n.spatial.units || ncol(dm) != n.spatial.units){
        stop(paste("Each distance matrix must be an n by n matrix. Spatial Units: ", n.spatial.units, ", dim(DM):", nrow(dm), ",", ncol(dm)))
     }
  }
  if (length(N) != n.spatial.units){
    stop(paste("N must be a vector of length equal to the number of spatial units (", n.spatial.units, ")", sep = ""))
  }

  if (length(spatial.factor) != nrow(x) || any(is.na(spatial.factor))){
    stop("spatial.factor must be the same length as the response, and can not contain missing values.")
  }
  if (var(table(spatial.factor)) != 0){
      stop("data must contain the same number of observations for each spatial unit, ordered by ascending time.")
  }
  n.time.points = NROW(I_star)/n.spatial.units
  if (n.time.points != floor(n.time.points)){
      # previous check should capture this case also... consider removing.
      stop("error when calculating number of time points - check your data.")  
  }
 
  if (verbose){
    cat("Building SEIR Model\n")
  }
  # Re-arrange N
  N = matrix(N, nrow = n.time.points, ncol = n.spatial.units, byrow=TRUE)

  # re-arrange I_star
  I_star_wide = matrix(0, nrow = n.time.points, ncol=n.spatial.units)
  for (i in 1:length(levels(spatial.factor))){
    I_star_wide[,i] = I_star[spatial.factor == levels(spatial.factor)[i]] 
  }
  I_star = I_star_wide

  # Calculate offset
  if (!missing(offset)) {
    if (length(offset) != n.time.points){ 
      stop(gettextf("number of offsets is %d, should equal %d (number of time points)", 
                    length(offset), n.time.points),
           domain = NA)
    }
  }
  else{
    cat("Assuming equally spaced case reports.\n")
    offset = rep(1, n.time.points)
  }
 
  # Calculate a reasonable starting value for the intercept
  priorBetaIntercept = log(mean(-log(1-(I_star/N))))
  
  dataModelInstance = buildDataModel(I_star, type="identity") 
  exposureModelInstance = buildExposureModel(X=x, 
                                             nTpt=n.time.points, 
                                             n.spatial.units,
                                             beta=c(priorBetaIntercept + rnorm(1), rnorm(ncol(x)-1)),
                                             betaPriorPrecision = rep(0.1, ncol(x)),
                                             betaPriorMean = rep(0, ncol(x)),
                                             offset = offset)
  reinfectionModelInstance = buildReinfectionModel("SEIR")
  samplingControlInstance = buildSamplingControl(iterationStride=1000) 
  distanceModelInstance = buildDistanceModel(distance.list) 
  transitionPriorsInstance = buildTransitionPriorsFromProbabilities(p_ei, p_ir, p_ei_ess, p_ir_ess) 
  
  N = matrix(N, ncol = n.spatial.units, nrow = NROW(I_star))
  E0 = I_star[1,]
  I0 = I_star[1,] + I_star[2,]
  S0 = N[1,] - E0 - I0
  
 
  initContainerInstance = buildInitialValueContainer(I_star, N, S0=S0,
                                                     I0=I0,E0=E0,
                                                     reinfection=FALSE,dataType="I_star") 
  modelObject = buildSEIRModel(chainFileName,dataModelInstance,exposureModelInstance,reinfectionModelInstance,distanceModelInstance,
                               transitionPriorsInstance, initContainerInstance, samplingControlInstance)
  modelObject$setRandomSeed(seed+1)
  modelObject$setTrace(0) 
  return(modelObject)
}

node.qSpatialSEIR = function(params){
    set.seed(seed + params[[1]])
    localModelObject <<- qSpatialSEIR(I_star ~ x,N,
                                      spatial.factor,
                                      distance.list,
                                      verbose, p_ei,
                                      p_ir, transition_ess,
                                      filenames[params[[1]]], seed + 10*params[[1]],
                                      offset=offset)
    for (i in 1:200)
    {
        localModelObject$simulate(10)
        localModelObject$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    for (i in 1:50)
    {
        localModelObject$simulate(100)
        localModelObject$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    localModelObject$parameterSamplingMode = 7
    localModelObject$compartmentSamplingMode = 17
    localModelObject$useDecorrelation = 10
    localModelObject$performHybridStep = 10
}

fit.qSpatialSEIR = function(formula, N, spatial.factor, distance.list=NA, verbose=TRUE, p_ei=NA, p_ir=NA, transition_ess=NA, seed=NA, n.cores=3, conv.criterion = 1.2, return.cluster=FALSE, data, offset)
{
    call <- match.call()
    if (verbose){
        cat("Processing user input.\n")
        cat("Call: ")
        print(call)
     cat("\n")
    }
    if (any(is.na(p_ei)) | any(is.na(p_ir))){
        warning("Not specifying meaningful prior values for the E-I or I-R transition probabilities can result in model fitting problems.")
        if(any(is.na(p_ei))){
            p_ei = 0.5
            p_ei_ess = 1
        }
        if(any(is.na(p_ir))){
            p_ir = 0.5
            p_ir_ess = 1
        }
    }
    else{
        if (all(is.na(transition_ess))){
            p_ei_ess = 1000
            p_ir_ess = 1000
        }
        else{
            p_ei_ess = transition_ess 
            p_ir_ess = transition_ess 
        }
    }
    if (all(is.na(seed))){
        seed = as.numeric(Sys.time())
        set.seed(seed)
    }
    if (missing(data))
    {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)  
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    cat("\n")
    I_star = model.response(mf, "numeric")
    mt <- attr(mf, "terms")

    if (is.empty.model(mt)){
        stop("Error, model formula was empty.")
    }  
    x <- model.matrix(mt, mf)

    # Process distance matrices and spatial.factor
    if (class(spatial.factor) != "factor"){
        spatial.factor = as.factor(spatial.factor)
    }

    if (all(is.na(distance.list)) || class(distance.list) != "list")
    {
        warning("Missing or invalid list of distance matrices. Assuming a single global spatial dependence term.")
        DM = (1-diag(length(levels(spatial.factor))))/(length(levels(spatial.factor)))
        distance.list = list(DM)
    }
    n.spatial.units = length(table(spatial.factor))
    for (dm in distance.list){
        if (class(dm) != "matrix" || nrow(dm) != n.spatial.units || ncol(dm) != n.spatial.units){
            stop(paste("Each distance matrix must be an n by n matrix. Spatial Units: ", n.spatial.units, ", dim(DM):", nrow(dm), ",", ncol(dm)))
        }
    }
    if (length(N) != n.spatial.units){
        stop(paste("N must be a vector of length equal to the number of spatial units (", n.spatial.units, ")", sep = ""))
    }

    if (length(spatial.factor) != nrow(x) || any(is.na(spatial.factor))){
        stop("spatial.factor must be the same length as the response, and can not contain missing values.")
    }
    if (var(table(spatial.factor)) != 0){
        stop("data must contain the same number of observations for each spatial unit, ordered by ascending time.")
    }
    n.time.points = NROW(I_star)/n.spatial.units
    if (n.time.points != floor(n.time.points)){
        # previous check should capture this case also... consider removing.
        stop("error when calculating number of time points - check your data.")  
    }

    if (verbose){
        cat("Building SEIR Model\n")
    }

    # Calculate offset
    if (!missing(offset)) {
        if (length(offset) != n.time.points){ 
            stop(gettextf("number of offsets is %d, should equal %d (number of time points)", 
                length(offset), n.time.points),
            domain = NA)
        }
    }
    else{
        cat("Assuming equally spaced case reports.\n")
        offset = rep(1, n.time.points)
    }

    # Build chains
    cl = makeCluster(n.cores)
    if (all(is.na(seed))){
        seed = as.numeric(System.time()) 
    }
    filenames = c("qSpSEIRchain1.csv", "qSpSEIRchain2.csv", "qSpSEIRchain3.csv")
    clusterExport(cl, c("I_star", "n.time.points", "n.spatial.units", 
                        "spatial.factor",
                        "distance.list", "x", "N", "verbose",
                        "p_ei", "p_ir", "transition_ess", "offset",
                        "seed", "filenames"), envir=environment())
    if (verbose){
        cat("Building chains.\n")
    }

    nodeParams = list(list(1), list(2), list(3))


    # build and burn in a model in each node
    parLapply(cl, nodeParams, node.qSpatialSEIR)
    conv=FALSE
    while(!conv){
        parLapply(cl, 1:3, function(x){localModelObject$simulate(25000)})
        conv = checkConvergence(filenames[1], filenames[2], filenames[3], 
                                maxVal = conv.criterion, verbose=verbose)
    }
    if (verbose){
        cat(paste("Chains converged, criterion: ", conv.criterion,"\n", sep = ""))
    }
    dat1=read.csv(filenames[1])
    dat2=read.csv(filenames[2])
    dat3=read.csv(filenames[3])
    if (!return.cluster){
        stopCluster(cl)
        return(list(mcmc.chain.1 = dat1,
                    mcmc.chain.2 = dat2,
                    mcmc.chain.3 = dat3))
    }
    else if (conv){
        return(list(mcmc.chain.1 = dat1,
                    mcmc.chain.2 = dat2,
                    mcmc.chain.3 = dat3,
                    cluster=cl))
    }
    else{
        cat("Sampling ended prematurely, killing cluster.\n")
        stopCluster(cl)
        return(list(mcmc.chain.1 = dat1,
                    mcmc.chain.2 = dat2,
                    mcmc.chain.3 = dat3))
    }

}
