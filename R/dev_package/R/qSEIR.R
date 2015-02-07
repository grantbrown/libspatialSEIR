qSEIR <- function(formula, N, verbose=TRUE, p_ei=NA, p_ir=NA, transition_ess=NA, chainFileName=NA, seed=NA, data, offset){
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
  if (length(N) != 1){
    warning("N was not of length 1, using only the first value.")
    N = N[1]
  }
  if (missing(data))
  {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)  
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  cat("\n")
  I_star = model.response(mf, "numeric")
  mt <- attr(mf, "terms")
  
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(I_star)){ 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    length(offset), NROW(y)),
           domain = NA)
    }
  }
  else{
    cat("Assuming equally spaced case reports.\n")
    offset = rep(1, NROW(I_star))
  }
  if (is.empty.model(mt)){
    stop("Error, model formula was empty.")
  }
  
  x <- model.matrix(mt, mf)
  
  if (verbose){
    cat("Building SEIR Model\n")
  }
  
  # Calculate a reasonable starting value for the intercept
  priorBetaIntercept = log(mean(-log(1-(I_star/N))))
  
  dataModelInstance = buildDataModel(I_star, type="identity") 
  exposureModelInstance = buildExposureModel(X=x, 
                                             nTpt=nrow(x), 
                                             1,
                                             beta=c(priorBetaIntercept, rep(0, ncol(x)-1)),
                                             betaPriorPrecision = rep(0.1, ncol(x)),
                                             betaPriorMean = rep(0, ncol(x)),
                                             offset = offset)
  reinfectionModelInstance = buildReinfectionModel("SEIR")
  samplingControlInstance = buildSamplingControl(iterationStride=1000) 
  distanceModelInstance = buildDistanceModel(list(matrix(0)))
  
  transitionPriorsInstance = buildTransitionPriorsFromProbabilities(p_ei, p_ir, p_ei_ess, p_ir_ess) 
  
  E0 = I_star[1]
  I0 = I_star[1]
  S0 = N - E0 - I0
  
  N = matrix(N, ncol = 1, nrow = NROW(I_star))
  
  initContainerInstance = buildInitialValueContainer(I_star, N, S0=S0,
                                                     I0=I0,E0=E0,
                                                     reinfection=FALSE,dataType="I_star") 
  modelObject = buildSEIRModel(chainFileName,dataModelInstance,exposureModelInstance,reinfectionModelInstance,distanceModelInstance,
                               transitionPriorsInstance, initContainerInstance, samplingControlInstance)
  modelObject$setRandomSeed(seed+1)
  modelObject$setTrace(0) 
  return(modelObject)
}

node.qSEIR = function(params){
    set.seed(seed + params[[1]])
    localModelObject <<- qSEIR(I_star ~x,N,
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

fit.qSEIR = function(formula, N, verbose=TRUE, p_ei=NA, p_ir=NA, transition_ess=NA, seed=NA, n.cores=3, conv.criterion = 1.2, data, offset){

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
    if (length(N) != 1){
        warning("N was not of length 1, using only the first value.")
        N = N[1]
    }
    if (missing(data))
    {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)  
    m <- match(c("formula", "data", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    cat("\n")
    I_star = model.response(mf, "numeric")
    mt <- attr(mf, "terms")

    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
    if (length(offset) != NROW(I_star)){ 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    length(offset), NROW(y)),
           domain = NA)
    }
    }
    else{
        cat("Assuming equally spaced case reports.\n")
        offset = rep(1, NROW(I_star))
    }
    if (is.empty.model(mt)){
    stop("Error, model formula was empty.")
    }

    x <- model.matrix(mt, mf)

    if (verbose){
        cat("Building SEIR Model\n")
    }

    # Build chains
    cl = makeCluster(n.cores)
    if (all(is.na(seed))){
        seed = as.numeric(System.time()) 
    }
    filenames = c("qSEIRchain1.csv", "qSEIRchain2.csv", "qSEIRchain3.csv")
    clusterExport(cl, c("I_star", "x", "N", "verbose",
                        "p_ei", "p_ir", "transition_ess", "offset",
                        "seed", "filenames"), envir=environment())
    if (verbose){
        cat("Building chains.\n")
    }

    nodeParams = list(list(1), list(2), list(3))


    # build and burn in a model in each node
    parLapply(cl, nodeParams, node.qSEIR)
    conv=FALSE
    while(!conv){
        parLapply(cl, 1:3, function(x){localModelObject$simulate(25000)})
        conv = checkConvergence(filenames[1], filenames[2], filenames[3], 
                                maxVal = conv.criterion, verbose=verbose)
    }
    if (verbose){
        cat(paste("Chains converged, criterion: ", conv.criterion,"\n", sep = ""))
    }
    dat1=read.csv(filenams[1])
    dat2=read.csv(filenams[2])
    dat3=read.csv(filenams[3])
    stopCluster(cl)
    return(list(mcmc.chain.1 = dat1,
                mcmc.chain.2 = dat2,
                mcmc.chain.3 = dat3))
}
