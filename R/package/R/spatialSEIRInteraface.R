spatialSEIRModel = function(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar, Estar, Istar, Rstar, offset, X_, Z_, X_pRS_, DistMat_, rho_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, p_ei_, p_ir_, N_, outFile, iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths, reinfectionMode)
{
    interface = new( spatialSEIRInterface )
    interface$buildSpatialSEIRInterface(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar, Estar, Istar, Rstar,offset, X_, Z_, X_pRS_, DistMat_, rho_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, p_ei_, p_ir_, N_, outFile, iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths, reinfectionMode)
    return(interface)
}
