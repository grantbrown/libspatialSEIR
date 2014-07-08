spatialSEIRModel = function(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar, Estar, Istar, Rstar, offset, X_, Z_, X_pRS_, DistMat_, rho_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, gamma_ei_, gamma_ir_, N_, outFile, iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths, reinfectionMode)
{
    interface = new( spatialSEIRInterface )
    err = interface$buildSpatialSEIRInterface(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar, Estar, Istar, Rstar,offset, X_, Z_, X_pRS_, DistMat_, rho_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, gamma_ei_, gamma_ir_, N_, outFile, iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths, reinfectionMode)
    if (err < 0){
        rm(interface)
        stop("Errors building spatialSEIRModel\n")
    }
    return(interface)
}
