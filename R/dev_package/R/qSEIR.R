qSEIR <- function(formula, N, verbose=TRUE, data, offset){
    # The following is based on the built in 'lm' function, 
    # and handles the formula and data provided by the user.
    call <- match.call()
    if (verbose){
        cat("Processing user input.\n")
        cat("Call: ")
        print(call)
        cat("\n")
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

    print(mf)
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
    
    if (verbose){
        cat("Building SEIR Model\n")
    }
    
    dataModel = buildDataModel(I_star, type="identity") 

}
